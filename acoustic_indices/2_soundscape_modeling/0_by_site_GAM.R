# Pre-modeling analysis of acoustic indices
# links to help:
# https://noamross.github.io/mgcv-esa-workshop/ and https://noamross.github.io/gams-in-r-course/
# 1) GAM basics https://gsp.humboldt.edu/olm/R/05_03_GAM.html
# 2) GAM coef simulations https://www.maths.ed.ac.uk/~swood34/mgcv/check-select.pdf
# 3) GAM book pdf in Mendeley"
#       - pg182 mgcv description of model fits
#       - pg328 describes some prediciton
# 4) GL + GAM general description: https://christophm.github.io/interpretable-ml-book/extend-lm.html#gam
# 5) GLM vs GAM opinion: https://stats.stackexchange.com/questions/380426/when-to-use-a-gam-vs-glm
# 6) some GAM implementation: https://m-clark.github.io/generalized-additive-models/application.html
# 7) posterior estimate with GAM: https://stats.stackexchange.com/questions/190348/can-i-use-bootstrapping-to-estimate-the-uncertainty-in-a-maximum-value-of-a-gam
# 8) https://kevintshoemaker.github.io/NRES-746/Generalized%20Additive%20Models%20(GAMs).pdf
# 9) Random effects post: https://stats.stackexchange.com/questions/552880/by-group-random-effect-gam
# 10) Hierarchical paper: https://peerj.com/articles/6876/
#     - should each group (ARU) have a smoother? 
#     - should each group smoother be the same wilggliness?
#     - will smoothers for each group be similar to one another?
#         - probably for all indices either G, GS, or I

# MODELING NOTES
# June 8: inital modeling
# June 16: scaled covariates, including no RandEff, including log(wavs), accounting for anamalous points in diagnostics
# NDSI: (June16) complete
# ACI: (June16) complete
# ADI: (June16) complete: resids plot has no vlaues at y = 0 (looks like a bivalve)
# AEI:
# BI: 
# H: 
# Hs: 
# Ht: 
# M: 
# NDSI_A: 
# NDSI_B: 
# R: 
# rugo: 
# sfm:
# zcr_mean:

# Outlier notes:
# 1) rugo/Int: found s2lam018_190506 & s2lam018_190412 with extreme interference (97%) shoulkd be removed. Skipping, electronic signal throughout recording

# Libraries used
# data libs
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)

# modeling libs
library(mgcv)
library(gratia)

# custom additional functions
source('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/code/paper1-AcousticIndices/utility_fxs.R')

###########################################
####### read in data and data prep ########
wd = '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/'
abgqi_df = fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI.csv')) %>%
  mutate(ARU = substr(site, 4,5),
         ARUdevice = substr(site, 4,8)) %>%
  select(-contains('var')) %>%
  select(site, wavs, contains('mean'), ARU, -Unidentified_mean) %>%
  rename(Anthropophony = Anthropophony_mean, Biophony = Biophony_mean, Geophony = Geophony_mean,
         Quiet = Quiet_mean, Interference = Interference_mean)
  
  
indices_df = fread(paste0(wd, 'acoustic_indices_aggregation/averages/site_avg_acoustic_indices.csv'))

# corrected ABGQI data (May 31 2022)
# abgqi_df = fread(paste0(wd, 'paired_ARUs/corrected_data/corrected_site_abgqi_03June2022.csv'))
# indices_df = fread(paste0(wd, 'paired_ARUs/corrected_data/corrected_site_acousticindices_03June2022.csv'))
# df = merge(indices_df, abgqi_df, by = 'site')


# Total minutes pre-analyses
sum(abgqi_df$wavs) # 741061

# Clean sites that are known problems
# 1. insufficient number of recordings
# 2. Sites known to have constant static error
# functioning sites always had > 2.5 days (60 hours) of recordings
abgqi_df = abgqi_df %>%
  filter(wavs >= (144 * 2.5)) # 44 sites

# sites with static
# s2lam050_190411, s2lam049_210501 : a lot of static but still has other signal
static_sites = paste("s2lam018_190517", # found in ACI analysis
                     "s2lam018_190412", "s2lam018_190506", # found in Rugo analysis
                     sep = "|")

abgqi_df = abgqi_df[!grepl(static_sites, abgqi_df$site),] # 1 site still in list - the site with highest %Interfernce 

# Clean sites with anamolous values
# 1. Sites with any value == 1
# only one site with any value == 1 (quiet)
mod_df = abgqi_df %>%
  filter(Quiet < 1)

# Number of minutes in model data
sum(mod_df$wavs) #728766


# visualize recordings against acoustic measures
temp_df = merge(x = indices_df, y = mod_df) 
temp_df %>%
  select(-site, -YYYY) %>%
  gather(variable, value, -wavs, -ARU) %>%
  ggplot(aes(x = log(wavs), y = value)) +
    geom_point(alpha = 0.5) +
    facet_wrap(~variable, scales = "free")

# scale covariates ABGQI
mod_df = mod_df %>%
  mutate(Anthropophony = min_max_norm(Anthropophony),
         Biophony      = min_max_norm(Biophony),
         Geophony      = min_max_norm(Geophony),
         Quiet         = min_max_norm(Quiet),
         Interference  = min_max_norm(Interference))

# list to store model slopes
model_slopes = list()

###########################################
################# NDSI ####################
i = "NDSI"
temp_y = data.frame(site = indices_df$site, NDSI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) %>%
   select(-c(site))

# visualize data
temp_df %>%
  gather(soundType, value, -NDSI, -ARU) %>%
  ggplot(aes(x = value, y = NDSI, colour = ARU)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "gam") +
    facet_wrap(~soundType, scales = "free")
hist(temp_df$NDSI)

# add transformed NDSI
temp_df$beta_NDSI = (temp_df$NDSI + 1)/2
hist(temp_df$beta_NDSI)

# Outliers
# Visual inspection and model diagnositcs - no gross outliers that cause alarm

# visualize pairs of corr
# temp_df %>%
#   GGally::ggpairs(aes(alpha = 0.05), progress = FALSE)

###### TESTING using: https://www.youtube.com/watch?v=Ukfvd8akfco
# start with full model under beta distribution
mod1 = gam(beta_NDSI ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = betar(),
           method = 'ML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # looks good
#concurvity(mod1, full = F)
draw(mod1, scales = 'fixed')

# drop geophony (p = 0.998)
mod2 = gam(beta_NDSI ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = betar(),
           method = 'ML')
summary(mod2)
gam.check(mod2)
AIC(mod1, mod2)

# drop logwavs (p = 0.181)
mod3 = gam(beta_NDSI ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = betar(),
           method = 'ML')
summary(mod3)
gam.check(mod3)
AIC(mod2, mod3)
draw(mod3, scales = "fixed", ) # nice visualization tool from gratia for partial effects
# smoothing parameters are related to variance components - assesses variation associated with random effects
# proportion of the variance attributed to RE's main effect
#variance_comp(mod4)

# QQ should follow 1:1 wihtin the reference band (do data follow assumptions of model well) some do not
# residuals should be approx gaussian 
vis.gam(mod3, theta = 65)
gam.vcomp(mod3)

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

# visualize predictions vs observed
pred_plot(data_df = temp_df, model_fit = mod4, index_name = "beta_NDSI")

###########################################
################# ACI #####################
i = "ACI"
temp_y = data.frame(site = indices_df$site, ACI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -ACI, -ARU, -wavs, -site) %>%
  ggplot(aes(x = value, y = ACI, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = ACI, fill = ARU)) +
  geom_histogram()

# add transformed
# bring ACI down to 0
temp_df$ACI_0 = temp_df$ACI - min(temp_df$ACI) + 1e-7

# both ARUs have similar shape and wiggliness -- if Randeffect is included use GS model type
ggplot(temp_df, aes(x = ACI_0, fill = ARU)) +
  geom_histogram()

# Extreme ACI values (high):
# 528 (s2lam050_190411)

# fit a gamma model based on positive real 
mod1 = gam(ACI_0 ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = Gamma(link = 'log'),
           method = 'ML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # heavy tail in QQ - primarily de to one observation

resids = data.frame(mod1$residuals) # most extreme obs is the only ACI = 0 
temp_df = temp_df[-381,]

# remodel without extreme outlier
remod1 = gam(ACI_0 ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = Gamma(link = 'log'),
           method = 'ML')
summary(remod1)
par(mfrow = c(2, 2))
gam.check(remod1, rep = 500) # good diagnostics now
AIC(mod1, remod1)

# remove logwavs (p = 0.6685)
mod2 = gam(ACI_0 ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = Gamma(link = 'log'),
           method = 'ML')
summary(mod2)
gam.check(mod2, rep = 500)
AIC(remod1, mod2)

# remove Quiet (p = 0.16)
mod3 = gam(ACI_0 ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = Gamma(link = 'log'),
           method = 'ML')
summary(mod3)
AIC(mod2, mod3)
gam.check(mod3, rep = 500)


draw(mod3, scales = "fixed") 
appraise(mod3, method = "simulate")

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

# visualize predictions vs observed
pred_plot(data_df = temp_df, model_fit = mod3, index_name = "ACI_0")

###########################################
################# ADI #####################
i = "ADI"
temp_y = data.frame(site = indices_df$site, ADI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -ADI, -ARU) %>%
  ggplot(aes(x = value, y = ADI, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = ADI, fill = ARU)) +
  geom_histogram()

# add transformed 
temp_df$ADInorm = min_max_norm(temp_df$ADI) # bring to 0-1 scale to allow for beta error distribution
ggplot(temp_df, aes(x = ADInorm, fill = ARU)) +
  geom_histogram()

# fit a model with gamma (real positive)
mod1 = gam(ADI ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'ML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1, rep = 500) # heavy QQ plot tail and possible residual skew

# try with betar family distribution
mod2 = gam(ADInorm ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = betar(),
           method = 'ML')
summary(mod2)
gam.check(mod2, rep = 500, ) # some odd fit issues in QQ and resids - gaussian appears better option

# check on the extreme residuals
resids = data.frame(residuals.gam(mod2))
plot(x = temp_df$ADInorm, y = resids$mod2.residuals)

temp_df = temp_df[-c(1033),] # lowest ADInorm value

# remodel without extreme outlier
remod2 = gam(ADInorm ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = betar(),
           method = 'ML')
summary(remod2)
gam.check(remod2, rep = 500) # some odd fit issues in resids plot but are still roughly normal. Two low QQ tail deviations


# remove logwavs (p = 0.879)
mod3 = gam(ADInorm ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = betar(),
           method = 'ML')
summary(mod3)
AIC(mod2, mod3)
gam.check(mod3, rep = 500) # get more deviation in QQ tails now

draw(mod3, scales = "fixed")

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

# visualize predictions vs observed
pred_plot(data_df = temp_df, model_fit = mod3, index_name = "ADInorm")

###########################################
################# AEI #####################
i = "AEI"
temp_y = data.frame(site = indices_df$site, AEI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         AEI = min_max_norm(AEI)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -AEI, -ARU) %>%
  ggplot(aes(x = value, y = AEI, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = AEI, fill = ARU)) +
  geom_histogram()

# add transformed 
temp_df$sqrtAEI = sqrt(temp_df$AEI) # approximates normal distribution
ggplot(temp_df, aes(x = sqrtAEI, fill = ARU)) +
  geom_histogram()

# fit a model
mod1 = gam(AEI ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # model overall appears okay aside from hard line in residuals

# add random effect
mod2 = gam(AEI ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod2)
AIC(mod1, mod2) #random effect is suggested
gam.check(mod2)

draw(mod2, scales = "fixed")

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

###########################################
################# BI ######################
i = "BI"
temp_y = data.frame(site = indices_df$site, BI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         BI = min_max_norm(BI)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -BI, -ARU) %>%
  ggplot(aes(x = value, y = BI, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)

# both ARUs have similar distribution - start with gaussian or gamma
ggplot(temp_df, aes(x = BI, fill = ARU)) +
  geom_histogram()

# fit a model
mod1 = gam(BI ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1)
par(mfrow = c(2, 2)) 
gam.check(mod1) # gaussian probably not working for these data based on QQ

# gamma
mod2 = gam(BI ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod2)
gam.check(mod2) # Looks like gaussian is more appropriate

# Add in Random Effect
mod3 = gam(BI ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           method = 'REML')
summary(mod3)
AIC(mod1, mod3) # mod3 preferred 
gam.check(mod3) # increase functions for quiet

# remove Anthro p > 0.05
mod4 = gam(BI ~ 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           method = 'REML')
summary(mod4)
AIC(mod3, mod4) # mod4 preferred 
gam.check(mod4)

# increase quiet basis fxs
mod5 = gam(BI ~ 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet, k = 25) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           method = 'REML')
summary(mod5)
AIC(mod4, mod5) # mod4 preferred 
gam.check(mod5)


draw(mod4, scales = "fixed")

# check for response scale compared to raw data 
mean(fitted(mod1)); mean(fitted(mod3)); mean(fitted(mod4)); mean(temp_df$BI) 

pred_plot(data_df = temp_df, model_fit = mod4, index_name = i)

# summarized slopes
ci = slope_summary(mod4)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

###########################################
################# H #######################
i = "H"
temp_y = data.frame(site = indices_df$site, H = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         H = min_max_norm(H)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -H, -ARU) %>%
  ggplot(aes(x = value, y = H, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = H, fill = ARU)) +
  geom_histogram()
summary(temp_df$H)


# fit a model
mod1 = gam(H ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # QQ plot has lower tail issues

# try with betar family distribution
mod2 = gam(H ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod2)
gam.check(mod2, rep = 500) # better diagnostics, still some influential points though in QQ lower end, low shoulders in residuals histogram

# add ARU Random effect
mod3 = gam(H ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod3)
AIC(mod2, mod3)
gam.check(mod3) # approximately same QQ, residuals have slight taper pattern, histogram looks better though
draw(mod3, scales = "fixed")


# try heavy tail distribution
mod4 = gam(H ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = scat(link = "identity"), # allows for heavier tail distribution
           method = 'REML')
summary(mod4)
gam.check(mod4) # approximately same QQ, residuals have slight taper pattern, histogram looks better though
draw(mod4, scales = "fixed")

# even though diagnostics look possibly better using scat the deviations of predictions do not folow the data as well as gaussian
pred_plot(data_df = temp_df, model_fit = mod4, index_name = i)
pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

###########################################
################# Hs ######################
i = "Hs"
temp_y = data.frame(site = indices_df$site, Hs = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         Hs = min_max_norm(Hs)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -Hs, -ARU) %>%
  ggplot(aes(x = value, y = Hs, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = Hs, fill = ARU)) +
  geom_histogram()
summary(temp_df$Hs)

# add transformed 
# NA

# fit a model
mod1 = gam(Hs ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # QQ plot has lower tail issues


# try with betar family distribution
mod2 = gam(Hs ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod2)
gam.check(mod2) # better diagnostics - QQ appears off based only on ~2 points, low shoulders in residuals histogram

# add ARU Random effect
mod3 = gam(Hs ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod3)
AIC(mod2, mod3)
gam.check(mod3) # slightly worse QQ patterns, residuals have possible taper pattern, histogram is okay
draw(mod3, scales = "fixed")

pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

###########################################
################# Ht ######################
i = "Ht"
temp_y = data.frame(site = indices_df$site, Ht = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         Ht = min_max_norm(Ht)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -Ht, -ARU) %>%
  ggplot(aes(x = value, y = Ht, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = Ht, fill = ARU)) +
  geom_histogram()
summary(temp_df$Ht)

# add transformed 
# 

# fit a model
mod1 = gam(Ht ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # QQ plot has lower tail issues, also possible heteroskedasticity 

# try with betar family distribution
mod2 = gam(Ht ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod2)
gam.check(mod2) # better diagnostics, still some influential points though in QQ lower end

# add ARU Random effect
mod3 = gam(Ht ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod3)
AIC(mod2, mod3)
gam.check(mod3) # slightly worse QQ patterns, residuals have possible taper pattern, histogram still has sharp peak
draw(mod3, scales="fixed")

# Model 3
pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

###########################################
################# M ######################
i = "M"
temp_y = data.frame(site = indices_df$site, M = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         M = min_max_norm(M))

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -M, -ARU) %>%
  ggplot(aes(x = value, y = M, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = M, fill = ARU)) +
  geom_histogram()
summary(temp_df$M)

# Extreme high values are in M
# 73 - no discernable patterns
# 808, 524 - both have extremely high amounts of interference

# fit a model
mod1 = gam(M ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # wrong family in model - QQ plot plot has heavy deviations, increasing variance in residuals

# try beta family
mod2 = gam(M ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = betar,
           method = 'REML')
summary(mod2)
gam.check(mod2) # still some heavy deviations in QQ plot, improved heterosk in resids plot

# try gamma family
mod3 = gam(M ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod3)
gam.check(mod3) # looks like gamma family is the way to go - Q plot tails are brought down with gamma

# add random effect
mod4 = gam(M ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod4)
AIC(mod3, mod4)
gam.check(mod4)

pred_plot(data_df = temp_df, model_fit = mod4, index_name = i)

# summarized slopes
ci = slope_summary(mod4)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

###########################################
################# NDSI_A ##################
i = "NDSI_A"
temp_y = data.frame(site = indices_df$site, NDSI_A = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         NDSI_A = min_max_norm(NDSI_A))

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -NDSI_A, -ARU) %>%
  ggplot(aes(x = value, y = NDSI_A, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = NDSI_A, fill = ARU)) +
  geom_histogram()
summary(temp_df$NDSI_A)

# No extreme values in NDSI_A

# fit a model
mod1 = gam(NDSI_A ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # family could be okay

# beta
mod2 = gam(NDSI_A ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = betar,
           method = 'REML')
summary(mod2)
gam.check(mod2) # more extreme deviation from QQ

# gaussian + Random effects
mod3 = gam(NDSI_A ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           method = 'REML')
summary(mod3)
gam.check(mod3) # marginally better QQ pattern and residuals 
AIC(mod1, mod3)


# geophony p >0.05
mod4 = gam(NDSI_A ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod4)
gam.check(mod4) # no increase in edf
AIC(mod3, mod4) # mod4

pred_plot(data_df = temp_df, model_fit = mod4, index_name = i)

# summarized slopes
ci = slope_summary(mod4)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

###########################################
################# NDSI_B ##################
i = "NDSI_B"
temp_y = data.frame(site = indices_df$site, NDSI_B = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         NDSI_B = min_max_norm(NDSI_B))

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -NDSI_B, -ARU) %>%
  ggplot(aes(x = value, y = NDSI_B, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = NDSI_B, fill = ARU)) +
  geom_histogram()
summary(temp_df$NDSI_B)

# No extreme values in NDSI_B

# fit a model
mod1 = gam(NDSI_B ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # diagnostics look good

# Gaussian + Random effects
mod2 = gam(NDSI_B ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod2)
gam.check(mod2) # equivalent QQ pattern and residuals
AIC(mod1, mod2) # Random effect recommended

# Geophony p > 0.05
mod3 = gam(NDSI_B ~ 
             s(Anthropophony, k = 25) + 
             s(Biophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod3)
gam.check(mod3) 
AIC(mod2, mod3) # Mod2 with geophony suggested

pred_plot(data_df = temp_df, model_fit = mod2, index_name = i)

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

###########################################
################# R #######################
i = "R"
temp_y = data.frame(site = indices_df$site, R = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         R = min_max_norm(R))

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -R, -ARU) %>%
  ggplot(aes(x = value, y = R, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = R, fill = ARU)) +
  geom_histogram()
summary(temp_df$R) # positive continuous

# A few high, extreme values in R
# 223 (s2lam021_200710), 202 (s2lam019_200710), 365 (s2lam036_200710) - no discernible patterns

# add transformed 
temp_df$logR = log(temp_df$R)
ggplot(temp_df, aes(x = logR, fill = ARU)) +
  geom_histogram()

# fit a model
mod1 = gam(R ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1) # TERRIBLE performance
par(mfrow = c(2, 2))
gam.check(mod1) # Terrible approximation 

# gamma
mod2 = gam(R ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod2)
gam.check(mod2) # better diagnostics with extreme tails

# drop biophony (p > 0.05)
mod3 = gam(R ~ 
             s(Anthropophony) +
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
AIC(mod2, mod3) # dropping bio is good
summary(mod3)
gam.check(mod3) # better diagnostics


# doesnt work with slope offsets (June 13)
pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

###########################################
################# rugo ####################
i = "rugo"
temp_y = data.frame(site = indices_df$site, rugo = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         rugo = min_max_norm(rugo))

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -rugo, -ARU) %>%
  ggplot(aes(x = value, y = rugo, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = rugo, fill = ARU)) +
  geom_histogram()
summary(temp_df$rugo) # positive continuous

# No longer two extremely high values (June13)
# IDed and removed both (s2lam018_190506) and (s2lam018_190412) both very high interference

# fit a basic model
mod1 = gam(rugo ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1) 
par(mfrow = c(2, 2))
gam.check(mod1) # irregular upward concavity tails in QQ and residuals plots


# fit a Gamma structure model
mod2 = gam(rugo ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod2) 
gam.check(mod2) # better error structure in QQ wiht a few tail issues and variance in residuals.

# Gamma + Random effects
mod3 = gam(rugo ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod3)
gam.check(mod3) # alright diagnostics. Upper QQ is now pulled upward though. Probably the points in the response ~ fitted that deviates from 1:1
AIC(mod2, mod3) # Random effect recommended

pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

###########################################
################# sfm #####################
i = "sfm"
temp_y = data.frame(site = indices_df$site, sfm = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         sfm = min_max_norm(sfm))

# visualize data - definitely different y-intercepts for ARU
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -sfm, -ARU) %>%
  ggplot(aes(x = value, y = sfm, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
# approximately normal with different mean/sd 
ggplot(temp_df, aes(x = sfm, fill = ARU)) +
  geom_histogram()
summary(temp_df$sfm) # positive continuous

# No extreme values

# fit a basic model
mod1 = gam(sfm ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1) 
par(mfrow = c(2, 2))
gam.check(mod1)

# Add in Random effect first
mod2 = gam(sfm ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
gam.check(mod2) # QQ plots is okay, residuals has obvious ARU pattern and decreasing variance for AM
summary(mod2)

# Rm Anthro (largest p > 0.05)
mod3 = gam(sfm ~ 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
gam.check(mod3) # same diagnositc issue as mod2
summary(mod3)
AIC(mod2, mod3) # not delta3

# Rm Bio (largest p > 0.05)
mod4 = gam(sfm ~ 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
gam.check(mod4) # same diagnositc issue as mod2
summary(mod4)
AIC(mod2, mod4) #deltaAIC < 3 and mod4 is more parsimonious

pred_plot(data_df = temp_df, model_fit = mod4, index_name = i)

# summarized slopes
ci = slope_summary(mod4)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)

###########################################
################# zcr_mean ################
i = "zcr_mean"
temp_y = data.frame(site = indices_df$site, zcr_mean = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         zcr_mean = min_max_norm(zcr_mean))

# visualize data - definitely different y-intercepts for ARU
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -zcr_mean, -ARU) %>%
  ggplot(aes(x = value, y = zcr_mean, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = zcr_mean, fill = ARU)) +
  geom_histogram()
summary(temp_df$zcr_mean) # positive continuous

# No extreme values

# fit a basic model
mod1 = gam(zcr_mean ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod1) 
par(mfrow = c(2, 2))
gam.check(mod1) # good diagnostics - try ARU RE

# Add in Random effect first
mod2 = gam(zcr_mean ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
gam.check(mod2) # QQ plot deviates a bit more than without RE but model is highly recommended over no RE
summary(mod2)
AIC(mod1, mod2) 

draw(mod2, scales = "fixed")
pred_plot(data_df = temp_df, model_fit = mod2, index_name = i)

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)


###########################################
######### Visualize effects ###############

# calculate all slope values
slope_df = do.call("rbind", model_slopes)

# use middle 99% of x values to plot slopes to avoid extreme tail behavior
abgqi_limits = temp_df %>%
  select(Anthropophony, Biophony, Geophony, Quiet, Interference) %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(lower = quantile(value, 0.005),
            upper = quantile(value, 0.995))

# filter slope df based on percentile values
slope_df_filtered = slope_df %>%
  filter(case_when(var == 'Anthropophony' ~ data >= abgqi_limits$lower[1] & data <= abgqi_limits$upper[1],
                   var == 'Biophony'      ~ data >= abgqi_limits$lower[2] & data <= abgqi_limits$upper[2],
                   var == 'Geophony'      ~ data >= abgqi_limits$lower[3] & data <= abgqi_limits$upper[3],
                   var == 'Interference'  ~ data >= abgqi_limits$lower[4] & data <= abgqi_limits$upper[4],
                   var == 'Quiet'         ~ data >= abgqi_limits$lower[5] & data <= abgqi_limits$upper[5]))
slope_df_filtered %>%
  group_by(var) %>%
  summarise(min = min(data, na.rm = TRUE),
            max = max(data, na.rm = TRUE))

# summarize slope trends
slope_cis = slope_df_filtered %>%
  select(-smooth) %>%
  group_by(index, var) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) %>%
  mutate(index = factor(index))

# summarize slope trends using CI of mean
# slope_cis = slope_df %>%
#   select(-smooth, -se, -crit, -lower, -upper) %>%
#   group_by(index, var) %>%
#   summarise(lower025 = quantile(derivative, 0.025),
#             upper975 = quantile(derivative, 0.975),
#             median = quantile(derivative, 0.50)) %>%
#   mutate(index = factor(index))

# simple 95% CI error plot of summary
avg_slope_error_plot = function(slope_avg_df, temp_index) {
  temp_min = floor(min(slope_avg_df$lower))
  temp_max = ceiling(max(slope_avg_df$upper))
  gg = slope_avg_df %>%
    group_by(var) %>%
    filter(var == temp_index) %>%
    #mutate(index = fct_reorder(index, derivative)) %>%
    ggplot(aes(x = index, y = derivative)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ylim(c(temp_min, temp_max)) + # min max in slope df
    ylab("Slope") +
    xlab("Acoustic Index") +
    coord_flip() +
    facet_wrap(~var) +
    theme_bw()
  return(gg)
}

ggAnthro = avg_slope_error_plot(slope_cis, 'Anthropophony')
ggBio = avg_slope_error_plot(slope_cis, 'Biophony')
ggGeo = avg_slope_error_plot(slope_cis, 'Geophony')
ggQuiet = avg_slope_error_plot(slope_cis, 'Quiet')
ggInt = avg_slope_error_plot(slope_cis, 'Interference')
allplots = ggarrange(ggAnthro, ggBio, ggGeo, ggQuiet, ggInt, 
          ncol = 3, nrow = 2,
          labels = c("A", "B", "C", "D", "E"))
annotate_figure(allplots, top = text_grob("95% CI for first derivative (slopes) of GAM partial effects.
Note: should be interpreted alongside partial effects plots.", size = 14, face = "bold"),
                bottom = text_grob("Values reflect middle 99% of covariate range based on erroneous tail behavior.",
                                   size = 10))          
          
# ggexport(bxp, dp, lp, bxp, filename = "test.pdf",
#          nrow = 2, ncol = 1)

# more complex color coded by x-value
# could code by LULC or ARU too
slope_df_filtered %>%
  ggplot(aes(x = index, y = derivative)) +
  geom_point(aes(colour = data), size = 2, alpha = 0.5) +
  geom_errorbar(data = slope_cis, 
                aes(x = index, ymin = lower, ymax = upper), 
                width = 0.5, colour = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  coord_flip() +
  ylab("Slope") +
  xlab("Acoustic Index") +
  labs(title = 'Confidence intervals for first derivative (slopes) of GAM partial effects.') +
  facet_wrap(~var) +
  theme_bw()




#####################################






































# PRED VISUAL
# make the prediction, add this and a column of standard errors to the prediction data.frame.
# COULD USE THIS VISUAL AND SPECIFY EACH K VALUE ACCORDING TO MODEL
# pred = data.frame(
#   Anthropophony = temp_df$Anthropophony,
#   Biophony = temp_df$Biophony,
#   Geophony = temp_df$Geophony,  
#   Quiet = temp_df$Quiet,
#   Interference = temp_df$Interference,
#   BI = temp_df$BI,
#   predict(mod4, type = "response", se = TRUE)) %>%
#   gather(variable, value, -fit, -se.fit, -BI)
# 
# (ggAnthro = pred %>%
#     filter(variable == "Biophony") %>%
#     ggplot(aes(x = value)) +
#     geom_point(aes(y = BI), size = 1, alpha = 0.5) +
#     geom_line(aes(y = fit), alpha = 0.3) + # fit values
#     geom_smooth(aes(y = fit), colour = "red", alpha = 0.3, method = "gam", formula = y ~ s(x, k = 5))) # smoothed fit


# PRED VISUAL
# make the prediction, add this and a column of standard errors to the prediction data.frame.
# mod_pred = cbind(temp_df,
#                  predict(mod4,
#                          se.fit=TRUE, 
#                          type="response"))
# ( t = mod_pred %>%
#   select(-c(ACI,log_ACI)) %>%
#   gather(soundType, value, -ACI_0, -ARU, -fit, -se.fit) %>%
#   ggplot(aes(x = value, y = ACI_0, colour = ARU)) +
#     geom_point(alpha = 0.5) +
#     facet_wrap(~soundType) +
#     geom_ribbon(aes(ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit), x = value),
#                 alpha = 0.3, linetype = 0) )



# PREDICTIONS
library(itsadug)
# pred = get_predictions(mod4, 
#                        cond = list(Anthropophony = seq(min(temp_df$Anthropophony, na.rm = TRUE),
#                                                             max(temp_df$Anthropophony, na.rm = TRUE),
#                                                             length.out = nrow(temp_df)), se = TRUE))
# n = nrow(temp_df)
# pred = get_predictions(mod4, 
#                        cond = list(Anthropophony = seq(0, 1, length.out = 101), se = TRUE))
# 
# pred = pred %>%
#   mutate(BI = 1)
# 
# temp_df %>%
#   ggplot(aes(x = Anthropophony, y = BI)) +
#   geom_point(alpha = 0.3) +
#   geom_ribbon(data = pred, alpha = 0.4, aes(ymin = fit - CI, ymax = fit + CI), show.legend = F, fill = 'forestgreen') +
#   geom_line(data = pred, aes(y = fit), show.legend = FALSE, color = 'forestgreen')
# 



# generate predictions using covariate values
# https://www.maths.ed.ac.uk/~swood34/mgcv/check-select.pdf
# y_hat = predict(mod, type = 'lpmatrix')
# 
# # simulate coefficients drawing from posterior
# library(MASS)
# eq = seq(min(boot_df[,1]), max(boot_df[,1]), length = 1000)
# br = mvrnorm(1000, coef(fit), vcov(fit))
# max.eq = rep(NA, 1000)
# for (i in 1:1000){ 
#   fv = y_hat %*% br[i,]
#   max.eq[i] = eq[fv==max(fv)]
# }
# 
# ci = quantile(max.eq[complete.cases(max.eq)], c(.025,.975))
# ci

# PREDICTION OPTION
y_hat = predict.gam(mod4, type = 'response')
preds = data.frame(y = temp_df[,1], y_hat)
ggplot(preds, aes(x = y_hat, y = y)) +
  geom_point() +
  geom_abline(slope=1, intercept=0)



# # bootstrap resampling w/ replacement
# set.seed(143)
# n_boot = 3
# replicates = lapply(1:n_boot, function(x) temp_df[sample(nrow(temp_df), nrow(temp_df), replace = TRUE),])

# 
# summary(fit)
# gam.check(fit)
# 
# # diagnostics
# plot.gam(fit, residuals = TRUE, # residuals
#          shade = TRUE, # 95% CI for the mean shape of effects
#          shade.col = "lightblue", seWithMean = TRUE, pages = 1, pch =1, cex = 1, all.terms = TRUE)
# plot(fit1, pages = 1, all.terms = TRUE) # partial effects plots
# 
# 
# # small p-values indicate residuals are not randomly distributed (i.e., not enough basis fxs)
# gam.check(fit1)
# # if p is significant, try increasing k: gam(y ~ s(x1, k = 12))
# 
# # values over ~0.8 for worst may mean issues with collinearity
# concurvity(fit, full = TRUE)
# concurvity(fit, full = FALSE) # use to id where the relationships are UNID and Int is highest
# 
# AIC(fit) # 1074.22
# AIC(fit1) # 1056.504
# 

