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
# 11) Deviance residuals: https://bookdown.org/ltupper/340f21_notes/deviance-and-residuals.html


# MODELING NOTES
# June 8: inital modeling
# June 16: scaled covariates, including no RandEff, including log(wavs), accounting for anamalous points in diagnostics
# NDSI: (June16) complete
# ACI: (June16) complete
# ADI: (June16) complete: resids plot has no vlaues at y = 0 (looks like a bivalve)
# AEI: (June17) complete 
# BI: (June17) complete
# H: (June17) complete
# Hs: (June17) complete
# Ht: (June17) complete
# M: (June17) complete
# NDSI_A: (June17) complete
# NDSI_B: (June17) complete
# R: (June20) complete
# rugo: (June20) complete
# sfm: (June20) complete
# zcr_mean: (June20) complete

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

# Year summary
abgqi_df$YY = substr(abgqi_df$site, 10, 11)
abgqi_df %>%
  group_by(YY) %>%
  summarise(n())

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
error_sites = read.csv('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/data/identified_problem_sites_witherrors.csv')
error_sites = paste0(error_sites$SiteID, collapse = '|')
abgqi_df_no_error = abgqi_df[!grepl(error_sites, abgqi_df$site),] # 7 sites dropped after filter above

# Number of minutes in model data
mod_df = abgqi_df_no_error
sum(mod_df$wavs) #726801


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

# save list of sites used in analyses
# sites = temp_df$site
# write.csv(data.frame(sites), 
#           paste0(wd, 'sites_used_in_GAMs_21July22.csv'), 
#           row.names = F)

# list to store model slopes
model_slopes = list()
final_models = list()

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

# drop geophony (p = 0.927)
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

# drop logwavs (p = 0.170)
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
final_models[[i]] = mod3

# visualize predictions vs observed
pred_plot(data_df = temp_df, model_fit = mod3, index_name = "beta_NDSI")

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
# s2lam050_190411

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
gam.check(mod1) # heavy tail in QQ - primarily due to one observation

resids = data.frame(mod1$residuals) # most extreme obs is the only ACI = 0 
temp_df_remod = temp_df[-381,]

# remodel without extreme outlier
remod1 = gam(ACI_0 ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = Gamma(link = 'log'),
           method = 'ML')
summary(remod1)
par(mfrow = c(2, 2))
gam.check(remod1, rep = 500) # good diagnostics now
AIC(mod1, remod1)

# remove logwavs (p = 0.668350)
mod2 = gam(ACI_0 ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = Gamma(link = 'log'),
           method = 'ML')
summary(mod2)
gam.check(mod2, rep = 500)
AIC(remod1, mod2)

# remove Quiet (p = 0.223)
mod3 = gam(ACI_0 ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
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
final_models[[i]] = mod3

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
gam.check(mod2, rep = 500, ) # some odd fit issues in QQ and resids

# check on the extreme residuals
resids = data.frame(residuals.gam(mod2))
plot(x = temp_df$ADInorm, y = resids$mod2.residuals)
temp_df_remod = temp_df[-c(1030),] # lowest ADInorm value

# remodel without extreme outlier
remod2 = gam(ADInorm ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = betar(),
           method = 'ML')
summary(remod2)
gam.check(remod2, rep = 500) # some odd fit issues in resids plot but are still roughly normal. Two low QQ tail deviations


# remove logwavs (p = 0.859)
mod3 = gam(ADInorm ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = betar(),
           method = 'ML')
summary(mod3)
AIC(mod2, mod3)
gam.check(mod3, rep = 500)

draw(mod3, scales = "fixed")

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod3

# visualize predictions vs observed
pred_plot(data_df = temp_df, model_fit = mod3, index_name = "ADInorm")

###########################################
################# AEI #####################
i = "AEI"
temp_y = data.frame(site = indices_df$site, AEI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -AEI, -ARU) %>%
  ggplot(aes(x = value, y = AEI, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = AEI, fill = ARU)) +
  geom_histogram()

# add transformed 
temp_df$sqrtAEI = sqrt(temp_df$AEI) # approximates normal distribution
ggplot(temp_df, aes(x = sqrtAEI, fill = ARU)) +
  geom_histogram()

# fit a model with beta [0,1]
mod1 = gam(AEI ~ 
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
gam.check(mod1) # model overall appears okay aside from possible heteroskadasticity, lower QQ tail heavy

# remove logwavs (p = 0.897)
mod2 = gam(AEI ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = betar(),
           method = 'ML')
summary(mod2)
AIC(mod1, mod2) 
gam.check(mod2)

draw(mod2, scales = "fixed")

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod2

###########################################
################# BI ######################
i = "BI"
temp_y = data.frame(site = indices_df$site, BI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -BI, -ARU) %>%
  ggplot(aes(x = value, y = BI, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")

# both ARUs have similar distribution - start with gaussian or gamma
ggplot(temp_df, aes(x = BI, fill = ARU)) +
  geom_histogram()
summary(temp_df$BI)

# fit a model bound to positive real line
mod1 = gam(BI ~ 
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
gam.check(mod1) # diagnostics look okay - stick with Gamma(log)

# remove anthro (p = 0.517)
mod2 = gam(BI ~ 
             ARU +
             logwavs +
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'ML')
summary(mod2)
gam.check(mod2) # Looks like gaussian is more appropriate

# remove logwqavs (p = 0.524)
mod3 = gam(BI ~ 
             ARU +
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'ML')
summary(mod3)
AIC(mod2, mod3) # mod3 preferred 
gam.check(mod3)

draw(mod3, scales = "fixed")

# check for response scale compared to raw data 
mean(fitted(mod1)); mean(fitted(mod2)); mean(fitted(mod3)); mean(temp_df$BI) 

pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod3

# visualize predictions vs observed
pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

# clean up models
rm(mod1); rm(mod2); rm(mod3); rm(mod4); rm(mod5); rm(remod1); rm(remod2)

###########################################
################# H #######################
i = "H"
temp_y = data.frame(site = indices_df$site, H = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -H, -ARU) %>%
  ggplot(aes(x = value, y = H, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = H, fill = ARU)) +
  geom_histogram()
summary(temp_df$H)


# fit a model bnound to 0-1
mod1 = gam(H ~ 
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
gam.check(mod1) # QQ plot has lower tail deviation, resids decrease in var but are normal overall

# remove logwavs (p = 0.859)
mod2 = gam(H ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = betar(),
           method = 'ML')
summary(mod2)
gam.check(mod2, rep = 500)
AIC(mod1, mod2)


# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod2

pred_plot(data_df = temp_df, model_fit = mod2, index_name = i)

# clean up models
rm(mod1); rm(mod1_); rm(mod2); rm(mod3); rm(mod3_); rm(mod4); rm(mod5);rm(remod1);rm(remod2)

###########################################
################# Hs ######################
i = "Hs"
temp_y = data.frame(site = indices_df$site, Hs = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -Hs, -ARU) %>%
  ggplot(aes(x = value, y = Hs, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = Hs, fill = ARU)) +
  geom_histogram()
summary(temp_df$Hs)

# add transformed 
# NA

# fit a model bound 0-1
mod1 = gam(Hs ~ 
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
gam.check(mod1) # QQ plot has lower tail issues but otherwise good


# drop logwavs (p=0.866)
mod2 = gam(Hs ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = betar(),
           method = 'ML')
summary(mod2)
gam.check(mod2) # better diagnostics - QQ appears off based only on ~2 points, low shoulders in residuals histogram


pred_plot(data_df = temp_df, model_fit = mod2, index_name = i)

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod2

# clean up models
rm(mod1); rm(mod2); rm(mod3); rm(mod4); rm(mod5);rm(remod1);rm(remod2)

###########################################
################# Ht ######################
i = "Ht"
temp_y = data.frame(site = indices_df$site, Ht = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -Ht, -ARU) %>%
  ggplot(aes(x = value, y = Ht, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scale = "free")
ggplot(temp_df, aes(x = Ht, fill = ARU)) +
  geom_histogram()
summary(temp_df$Ht)

# add transformed 
# NA

# fit a model
mod1 = gam(Ht ~ 
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
gam.check(mod1) # good diags - maybe pronounced residual shoulders too

# drop logwavs (p = 0.332)
mod2 = gam(Ht ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = betar(),
           method = 'ML')
summary(mod2)
gam.check(mod2) # residuals slightly shifted to left skew
AIC(mod1, mod2)


# Model 2
pred_plot(data_df = temp_df_remod, model_fit = mod2, index_name = i)

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod2

rm(mod1); rm(mod2)

###########################################
################# M ######################
i = "M"
temp_y = data.frame(site = indices_df$site, M = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs))

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -M, -ARU) %>%
  ggplot(aes(x = value, y = M, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = M, fill = ARU)) +
  geom_histogram()
summary(temp_df$M)

# Extreme high values are in M
# 73 - no discernable patterns
# 808, 524 - both have extremely high amounts of interference
# transform
temp_df$normM = min_max_norm(temp_df$M)
temp_df$logM = log(temp_df$M)
ggplot(temp_df, aes(x = normM, fill = ARU)) +
  geom_histogram()
ggplot(temp_df, aes(x = logM, fill = ARU)) +
  geom_histogram()

# fit a model
mod1 = gam(M ~ 
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
gam.check(mod1) # Issues in diagnostics. heavy QQ tails and right skewed resids

# try normalized data
mod2 = gam(normM ~ 
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
gam.check(mod2) # better diagnostics except for one extreme point

# check on the extreme residuals
resids = data.frame(residuals.gam(mod2))
plot(x = temp_df$normM, y = resids$residuals.gam.mod2.)
temp_df_remod = temp_df[-c(73),] # highest norm M value

# model without largest outlier
remod2 = gam(normM ~ 
               ARU +
               logwavs +
               s(Anthropophony, k = 5) + 
               s(Biophony, k = 5) + 
               s(Geophony, k = 5) + 
               s(Quiet, k = 5) + 
               s(Interference, k = 5),
             data = temp_df_remod,
             family = betar(),
             method = 'ML')
summary(remod2)
par(mfrow = c(2, 2))
gam.check(remod2) # still some issues in diags

# try gaussian with log(M) with all data
mod3 = gam(logM ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           family = gaussian,
           method = 'ML')
summary(mod3)
gam.check(mod3) # best diags and look okay overall


# remove logwavs (p = 0.748)
mod4 = gam(logM ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = gaussian,
           method = 'ML')
summary(mod4)
gam.check(mod4)

# remove quiet (p = 0.113)
mod5 = gam(logM ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = gaussian,
           method = 'ML')
summary(mod5)
gam.check(mod5)


pred_plot(data_df = temp_df_remod, model_fit = mod5, index_name = "logM")

# summarized slopes
ci = slope_summary(mod5)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod5

rm(mod1); rm(mod2); rm(mod3); rm(mod4); rm(mod5); rm(remod2);
###########################################
################# NDSI_A ##################
i = "NDSI_A"
temp_y = data.frame(site = indices_df$site, NDSI_A = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs))

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -NDSI_A, -ARU) %>%
  ggplot(aes(x = value, y = NDSI_A, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = NDSI_A, fill = ARU)) +
  geom_histogram()
summary(temp_df$NDSI_A)

# No extreme values in NDSI_A

# fit a model
mod1 = gam(NDSI_A ~ 
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
gam.check(mod1) # good diags

# drop geophony (p = 0.613)
mod2 = gam(NDSI_A ~ 
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
gam.check(mod2) # more extreme deviation from QQ

# drop logwavs (p = 0.407)
mod3 = gam(NDSI_A ~ 
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
AIC(mod1, mod3)

pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod3

# clean up models
rm(mod1); rm(mod2); rm(mod3); 

###########################################
################# NDSI_B ##################
i = "NDSI_B"
temp_y = data.frame(site = indices_df$site, NDSI_B = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs))

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -NDSI_B, -ARU) %>%
  ggplot(aes(x = value, y = NDSI_B, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = NDSI_B, fill = ARU)) +
  geom_histogram()
summary(temp_df$NDSI_B)

# No extreme values in NDSI_B

# transforms
temp_df$normNDSI_B = min_max_norm(temp_df$NDSI_B)
ggplot(temp_df, aes(x = normNDSI_B, fill = ARU)) +
  geom_histogram()

# fit a model
mod1 = gam(NDSI_B ~ 
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
gam.check(mod1) # diagnostics are okay - QQ has deviating tails and resids are normal but decrease var

# Try norm data with beta
mod2 = gam(normNDSI_B ~ 
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
gam.check(mod2) # better diags with one larger deviance residual


# check on the extreme residuals
resids = data.frame(residuals.gam(mod2))
plot(x = temp_df$normNDSI_B, y = resids$residuals.gam.mod2.)
temp_df_remod = temp_df[-c(618),] # highest norm value

# remod mod2
remod2 = gam(normNDSI_B ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = betar(),
           method = 'ML')
summary(remod2)
par(mfrow = c(2, 2))
gam.check(remod2) # diags look good

# drop logwavs (p = 0.952)
mod3 = gam(normNDSI_B ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) +
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = betar(),
           method = 'ML')
summary(mod3)
gam.check(mod3) 
AIC(remod2, mod3)

# drop geophony (p = 0.210322)
mod4 = gam(normNDSI_B ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = betar(),
           method = 'ML')
summary(mod4)
gam.check(mod4) 
AIC(mod3, mod4)

pred_plot(data_df = temp_df, model_fit = mod4, index_name = "normNDSI_B")

# summarized slopes
ci = slope_summary(mod4)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod4

# clean up models
rm(mod1); rm(mod2); rm(mod3); rm(mod4); rm(remod2)

###########################################
################# R #######################
i = "R"
temp_y = data.frame(site = indices_df$site, R = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs))

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -R, -ARU) %>%
  ggplot(aes(x = value, y = R, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = R, fill = ARU)) +
  geom_histogram()
summary(temp_df$R) # positive continuous

# A few high, extreme values in R
# 223 (s2lam021_200710), 202 (s2lam019_200710), 365 (s2lam036_200710) - no discernible patterns

# add transformed 
temp_df$logR = log(temp_df$R)
ggplot(temp_df, aes(x = logR, fill = ARU)) +
  geom_histogram()

# fit a model for positive continuous values
mod1 = gam(R ~ 
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
gam.check(mod1) # OKay approximation, some QQ tail deviations

# try normal with log transform

mod2 = gam(logR ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           method = 'ML')
summary(mod2)
gam.check(mod2) # marginally better diagnostics

# drop geophony (p = 0.557679)
mod3 = gam(logR ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           method = 'ML')
AIC(mod2, mod3)
summary(mod3)
gam.check(mod3)

# drop logwavs (p = 0.15928)
mod4 = gam(logR ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           method = 'ML')
AIC(mod3, mod4)
summary(mod4)
gam.check(mod4)


pred_plot(data_df = temp_df, model_fit = mod4, index_name = 'logR')

# summarized slopes
ci = slope_summary(mod4)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod4

rm(mod1); rm(mod2); rm(mod3); rm(mod4);
###########################################
################# rugo ####################
i = "rugo"
temp_y = data.frame(site = indices_df$site, rugo = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs))

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -rugo, -ARU) %>%
  ggplot(aes(x = value, y = rugo, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = 'free')
ggplot(temp_df, aes(x = rugo, fill = ARU)) +
  geom_histogram()
summary(temp_df$rugo) # positive continuous

# No longer two extremely high values (June13)
# IDed and removed both (s2lam018_190506) and (s2lam018_190412) both very high interference

# fit a model with beta (positive [0,1])
mod1 = gam(rugo ~ 
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
gam.check(mod1) # decent fit with a possible outlier


# check on the extreme residuals
resids = data.frame(residuals.gam(mod1))
plot(x = temp_df$rugo, y = resids$residuals.gam.mod1.)
temp_df_remod = temp_df[-c(995),] # no systematic reason for higher residual

# refit model1
remod1 = gam(rugo ~ 
               ARU +
               logwavs +
               s(Anthropophony, k = 5) + 
               s(Biophony, k = 5) + 
               s(Geophony, k = 5) + 
               s(Quiet, k = 5) + 
               s(Interference, k = 5),
             data = temp_df_remod,
             family = betar(),
             method = 'ML')
summary(remod1) 
par(mfrow = c(2, 2))
gam.check(remod1) # structure is improved by dropping the outlier

# drop logwavs (p = 0.0589)
mod2 = gam(rugo ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = betar(),
           method = 'ML')
summary(mod2) 
AIC(remod1, mod2)
gam.check(mod2)

pred_plot(data_df = temp_df, model_fit = mod2, index_name = i)

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod2

rm(mod1); rm(mod2); rm(remod1);

###########################################
################# sfm #####################
i = "sfm"
temp_y = data.frame(site = indices_df$site, sfm = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs))

# visualize data - definitely different y-intercepts for ARU
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -sfm, -ARU) %>%
  ggplot(aes(x = value, y = sfm, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")

# approximately normal 0-1 bound or postive continuous
ggplot(temp_df, aes(x = sfm, fill = ARU)) +
  geom_histogram()
summary(temp_df$sfm)

# No extreme values

# fit a model for positive continuous values
mod1 = gam(sfm ~ 
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
gam.check(mod1) # diagnostics are okay, a couple of possible outliers in QQ

# check on the extreme residuals
resids = data.frame(residuals.gam(mod1))
plot(x = temp_df$sfm, y = resids$residuals.gam.mod1.)
temp_df_remod = temp_df[-c(1085),] # lowest value

# refit without outlier
remod1 = gam(sfm ~ 
               ARU +
               logwavs +
               s(Anthropophony, k = 5) + 
               s(Biophony, k = 5) + 
               s(Geophony, k = 5) + 
               s(Quiet, k = 5) + 
               s(Interference, k = 5),
             data = temp_df_remod,
             family = Gamma(link = "log"),
             method = 'ML')
summary(remod1) 
par(mfrow = c(2, 2))
gam.check(remod1) # better diagnostics with smaller tail in QQ

# drop logwavs (p = 0.458)
mod2 = gam(sfm ~ 
               ARU +
               s(Anthropophony, k = 5) + 
               s(Biophony, k = 5) + 
               s(Geophony, k = 5) + 
               s(Quiet, k = 5) + 
               s(Interference, k = 5),
             data = temp_df_remod,
             family = Gamma(link = "log"),
             method = 'ML')
summary(mod2) 
gam.check(mod2)

# drop geophony (p = 0.11718)
mod3 = gam(sfm ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = Gamma(link = "log"),
           method = 'ML')
summary(mod3) 
gam.check(mod3)


pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod3

###########################################
################# zcr_mean ################
i = "zcr_mean"
temp_y = data.frame(site = indices_df$site, zcr_mean = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU),
         logwavs = log(wavs))

# visualize data - definitely different y-intercepts for ARU
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -zcr_mean, -ARU) %>%
  ggplot(aes(x = value, y = zcr_mean, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = zcr_mean, fill = ARU)) +
  geom_histogram()
summary(temp_df$zcr_mean) # positive continuous or [0-1]

temp_df$normZCR = min_max_norm(temp_df$zcr_mean)

# No extreme values

# fit a postive model
mod1 = gam(zcr_mean ~ 
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
gam.check(mod1) # okay diagnostics - try betar

# beta
mod2 = gam(normZCR ~ 
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
gam.check(mod2) # had to tell if diags are better - one extreme value

# check on the extreme residuals
resids = data.frame(residuals.gam(mod2))
plot(x = temp_df$normZCR, y = resids$residuals.gam.mod2.)
temp_df_remod = temp_df[-c(726),] # highest norm value

# refit mod2
remod2 = gam(normZCR ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df_remod,
           family = betar(),
           method = 'ML')
summary(remod2) 
par(mfrow = c(2, 2))
gam.check(remod2) # better assumption and marginally better diags with beta

# drop logwavs (p = 0.387)
mod3 = gam(normZCR ~ 
               ARU +
               s(Anthropophony, k = 5) + 
               s(Biophony, k = 5) + 
               s(Geophony, k = 5) + 
               s(Quiet, k = 5) + 
               s(Interference, k = 5),
             data = temp_df_remod,
             family = betar(),
             method = 'ML')
summary(mod3) 
gam.check(mod3)


draw(mod3, scales = "fixed")
pred_plot(data_df = temp_df, model_fit = mod3, index_name = "normZCR")

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod3

###########################################
######### SAVE OBJECTS ####################
saveRDS(final_models, file = paste0(wd, 'modeling/acoustic_indices/LM_model_objects/gam_model_objects_20June2022.RData'))
saveRDS(model_slopes, file = paste0(wd, 'modeling/acoustic_indices/LM_model_objects/gam_model_slopes_20June2022.RData'))

######### Visualize effects ###############



   
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

