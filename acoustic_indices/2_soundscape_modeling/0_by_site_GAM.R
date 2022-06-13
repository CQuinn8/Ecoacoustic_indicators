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
# NDSI: (June8) complete
# ACI: (June9) complete
# ADI: (June9) complete: check if the ARU re is okay when I include in Geophony
# AEI: (June9) complete
# BI: (June9) complete
# H: (June9) complete
# Hs: (June9) complete
# Ht: (June10) complete - if ADI re is okay try random slopes for geophony and Quiet
# M: (June10) complete
# NDSI_A: (June10) complete
# NDSI_B: (June10) complete - if ADI re is okay try random slopes for anthro and Quiet
# R: (June13) complete - low performance - if ADI re is okay try random slopes for anthro
# rugo: (June13) complete
# sfm: (June13) complete - strange residuals behavior with obvious grouping in linear predictor from LG ~ AM
# zcr_mean: (June13) complete

# Outlier notes:
# 1) rugo/Int: found s2lam018_190506 & s2lam018_190412 with extreme interference (97%) shoulkd be removed. Skipping, electronic signal throughout recording

# TODO:
# finalize model visual tool
# create results table with: family/link, R2, ABGQI trends

library(data.table)
library(dplyr)
library(caret)
library(ggplot2)
library(mgcv)
library(tidyverse)
library(outliers)
library(gratia)
source('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/code/paper1-AcousticIndices/utility_fxs.R')

wd = '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/'
abgqi_df = fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI.csv')) %>%
  select(-contains('var')) %>%
  select(site, wavs, contains('mean'), -Unidentified_mean) %>%
  mutate(ARU = substr(site, 4,5)) %>%
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
# 2. Zero averages -> median value
# only one site with any value == 1 (quiet)
mod_df = abgqi_df %>%
  filter(Quiet < 1)

# account for zero values
mod_df[mod_df == 0] = NA

# Which rows contain zeroes
zeros_df = mod_df[!complete.cases(mod_df),] # n = 20 all are quiet

# calculate medians
medians = mod_df %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)))

# input medians for 0 (now NA) values
mod_df = mod_df %>% 
  mutate(across(where(is.numeric), .fns = ~ifelse(is.na(.x), median(.x, na.rm=TRUE), .x)))

# Number of minutes in model data
sum(mod_df$wavs) #728766

###########################################
################# NDSI ####################
i = "NDSI"
temp_y = data.frame(site = indices_df$site, NDSI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU)) %>%
   select(-c(site, wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -NDSI, -ARU) %>%
  ggplot(aes(x = value, y = NDSI, colour = ARU)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "gam") +
    facet_wrap(~soundType)
hist(temp_df$NDSI)

# add transformed NDSI
temp_df$posNDSI = temp_df$NDSI + 1
temp_df$beta_NDSI = (temp_df$NDSI + 1)/2
temp_df$sq_beta_NDSI = (temp_df$beta_NDSI)^2
hist(temp_df$beta_NDSI)

# Outliers
# Visual inspection and model diagnositcs - no gross outliers that cause alarm

# visualize pairs of corr
# temp_df %>%
#   GGally::ggpairs(aes(alpha = 0.05), progress = FALSE)

###### TESTING using: https://www.youtube.com/watch?v=Ukfvd8akfco
mod1 = gam(NDSI ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
            data = temp_df)
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # looks okay - slight tail extreme values in QQ and maybe some nonconstatn variance

# model using beta transformed NDSI - bounded 0-1
mod2 = gam(beta_NDSI ~ 
            s(Anthropophony) + 
            s(Biophony) + 
            s(Geophony) + 
            s(Quiet) + 
            s(Interference),
          data = temp_df,
          family = betar(),
          method = 'REML')
summary(mod2) # model summary
gam.check(mod2) # diagnostics not drastically better than gaussian but bound on 0-1 assumptionis included here

# add ARU random intercept effect
mod3 = gam(beta_NDSI ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're'), # random intercept
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod3) # model summary
AIC(mod2, mod3)
gam.check(mod3)


# increase basis functions for Anthro and Quiet
mod4 = gam(beta_NDSI ~ 
             s(Anthropophony, k = 25) + 
             s(Biophony) + 
             s(Geophony) +
             s(Quiet, k = 25) + 
             s(Interference) +
             s(ARU, bs = 're'), # random intercept
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod4)
AIC(mod3, mod4)
AIC(mod2, mod4) # delta ARU
gam.check(mod4)

# model 4 is best option
draw(mod4, scales = "fixed") # nice visualization tool from gratia for partial effects

# smoothing parameters are related to variance components - assesses variation associated with random effects
# proportion of the variance attributed to RE's main effect
#variance_comp(mod4)

# QQ should follow 1:1 wihtin the reference band (do data follow assumptions of model well) some do not
# residuals should be approx gaussian 
vis.gam(mod4, theta = 65)
gam.vcomp(mod4)


###########################################
################# ACI #####################
i = "ACI"
temp_y = data.frame(site = indices_df$site, ACI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU))

# visualize data
temp_df %>%
  gather(soundType, value, -ACI, -ARU) %>%
  ggplot(aes(x = value, y = ACI, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = ACI, fill = ARU)) +
  geom_histogram()

# add transformed 
# bring ACI down to 0
temp_df$ACI_0 = temp_df$ACI - min(temp_df$ACI) + 1e-7
temp_df$log_ACI = log(temp_df$ACI)
temp_df$sqrt_ACI = sqrt(temp_df$ACI)

# both ARUs have similar shape and wiggliness -- if Randeffect is included use GS model type
ggplot(temp_df, aes(x = ACI_0, fill = ARU)) +
  geom_histogram()

# Extreme ACI values (high):
# 528 (s2lam050_190411)

# fit a normal model
mod1 = gam(ACI ~ 
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
gam.check(mod1) # Extreme upper QQ plot behavior, increasing variance in residuals too

# change family to gamma (log link) and use the 0 shifted ACI values
mod2 = gam(ACI_0 ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod2)
gam.check(mod2) # assumptions look much better - still one outlier that is pulling QQ and resids

# try adding random effect
mod3 = gam(ACI_0 ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod3)
AIC(mod2, mod3) # prefer mod3 with RE
gam.check(mod3, rep = 500)

# remove Quiet (p > 0.05)
mod4 = gam(ACI_0 ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod4)
AIC(mod3, mod4) # prefer mod4: parsimony
gam.check(mod4, rep = 500)

draw(mod4, scales = "fixed") 
appraise(mod4, method = "simulate") # only works with gaussian data

###########################################
################# ADI #####################
i = "ADI"
temp_y = data.frame(site = indices_df$site, ADI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU)) %>%
  select(-c(site,wavs))

# visualize data
temp_df %>%
  gather(soundType, value, -ADI, -ARU) %>%
  ggplot(aes(x = value, y = ADI, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType)
ggplot(temp_df, aes(x = ADI, fill = ARU)) +
  geom_histogram()

# add transformed 
temp_df$ADInorm = min_max_norm(temp_df$ADI, min_x = min(temp_df$ADI), max_x = max(temp_df$ADI)) # bring to 0-1 scale
ggplot(temp_df, aes(x = ADInorm, fill = ARU)) +
  geom_histogram()

# no gross outliers

# fit a model
mod1 = gam(ADI ~ 
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
gam.check(mod1) # model on raw scale has non-constant variance of residuals

# try scaled ADI 0-1 
mod2 = gam(ADInorm ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod2)
gam.check(mod2, rep = 500)

# try with betar family distribution
mod3 = gam(ADInorm ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod3)
gam.check(mod3, rep = 500) # better diagnostics, still some influential points though in QQ lower end

# add ARU Random effect
mod4 = gam(ADInorm ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're'),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod4)
AIC(mod3, mod4) # suggests RandEff model
gam.check(mod4) # get more deviation in QQ tails now

# k is too low for geophony - higher basis is just wiggle fitting. Try random slopes for geophony first
mod5 = gam(ADInorm ~ 
             ARU +
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony, by = ARU) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're'),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod5)
AIC(mod4, mod5) # suggests random slopes for Geophony
gam.check(mod5) # basis increase recommended
draw(mod5, scales = "fixed")


# check for response scale compared to raw data 
mean(fitted(mod1));mean(temp_df$ADI) # approximates exactly
mean(fitted(mod5)); mean(fitted(mod4)); mean(temp_df$ADInorm) # is within a reasonable margin ( - 0.08)


###########################################
################# AEI #####################
i = "AEI"
temp_y = data.frame(site = indices_df$site, AEI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU)) %>%
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
gam.check(mod1)
# model overall appears okay aside from hard line in residuals

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

###########################################
################# BI ######################
i = "BI"
temp_y = data.frame(site = indices_df$site, BI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU)) %>%
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
summary(temp_df$BI)

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
gam.check(mod2) # QQ more appropriate - still one low value and ~5 upper values not approximated. all other diags look good

# Add in Random Effect
mod3 = gam(BI ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod3)
AIC(mod2, mod3) # mod3 preferred 
gam.check(mod3) # increase functions for quiet

# increase basis fxs for Quiet
mod4 = gam(BI ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet, k = 25) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod4)
AIC(mod3, mod4) # mod3 preferred 
gam.check(mod4)
draw(mod4, scales = "fixed")

# check for response scale compared to raw data 
mean(fitted(mod1)); mean(fitted(mod3)); mean(fitted(mod4)); mean(temp_df$BI) 

pred_plot(data_df = temp_df, model_fit = mod4, index_name = i)

###########################################
################# H #######################
i = "H"
temp_y = data.frame(site = indices_df$site, H = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU)) %>%
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

# add transformed 
# NA

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

pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

###########################################
################# Hs ######################
i = "Hs"
temp_y = data.frame(site = indices_df$site, Hs = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU)) %>%
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
gam.check(mod2) # better diagnostics, low shoulders in residuals histogram

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

###########################################
################# Ht ######################
i = "Ht"
temp_y = data.frame(site = indices_df$site, Ht = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU)) %>%
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
# NA - has approximate beta distribution (high alpha low beta)

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

# increase basis functions for geophony and quiet
mod4 = gam(Ht ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony, k = 25) + 
             s(Quiet, k = 25) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod4)
AIC(mod3, mod4)
gam.check(mod4) # fixed k issues

# Model 3 based on increased basis not changing model
draw(mod3, scales = "fixed")
pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

###########################################
################# M ######################
i = "M"
temp_y = data.frame(site = indices_df$site, M = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU))

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

# remove Quiet (p > 0.05)
mod5 = gam(M ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod5)
AIC(mod4, mod5)
gam.check(mod5)

pred_plot(data_df = temp_df, model_fit = mod5, index_name = i)


###########################################
################# NDSI_A ##################
i = "NDSI_A"
temp_y = data.frame(site = indices_df$site, NDSI_A = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU))

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
gam.check(mod1) # family could be okay - based on 0-1 bounding try beta

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
gam.check(mod2) # margnially better QQ pattern and residuals

# beta + Random effects
mod3 = gam(NDSI_A ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = betar,
           method = 'REML')
summary(mod3)
gam.check(mod3) # marginally better QQ pattern and residuals 
AIC(mod2, mod3)

# need to increase basis fxs for quiet
mod4 = gam(NDSI_A ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet, k = 25) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod4)
gam.check(mod4) # no increase in edf
AIC(mod3, mod4) # mod3 is preferred

pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

###########################################
################# NDSI_B ##################
i = "NDSI_B"
temp_y = data.frame(site = indices_df$site, NDSI_B = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU))

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
gam.check(mod1) # diagnostics look good. Based on positive distribution try gamma error

# gamma
mod2 = gam(NDSI_B ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod2)
gam.check(mod2) # worse diagnostics - stick with gaussian

# Gaussian + Random effects
mod3 = gam(NDSI_B ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod3)
gam.check(mod3) # equivalent QQ pattern and residuals
AIC(mod1, mod3) # Random effect recommended

# need to increase basis fxs for quiet
mod4 = gam(NDSI_B ~ 
             s(Anthropophony, k = 25) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet, k = 25) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod4)
gam.check(mod4) # no increase in edf
AIC(mod3, mod4) # no significant difference in models - mod3 can be used via parsimony

pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

###########################################
################# R #######################
i = "R"
temp_y = data.frame(site = indices_df$site, R = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU))

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
temp_df$Rnorm = min_max_norm(temp_df$R, min_x = min(temp_df$R), max_x = max(temp_df$R)) # bring to 0-1 scale
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
gam.check(mod2) # better diagnostics but only marginally - still high tails


# beta with normalized values
mod3 = gam(Rnorm ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod3)
gam.check(mod3) 


# Gaussian + Random effects
mod3 = gam(R ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
summary(mod3)
gam.check(mod3) # equivalent QQ pattern and residuals
AIC(mod1, mod3) # Random effect recommended

# transform R and use haby tailed dist version of normal
mod4 = gam(logR ~ 
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
gam.check(mod4) # better diagnostics

# drop geophony (p > 0.05)
mod5 = gam(logR ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = scat(link = "identity"),
           method = 'REML')
AIC(mod4, mod5)
summary(mod5)
gam.check(mod5) # better diagnostics

# Anthro Lg and AM have differnt behavior
mod6 = gam(logR ~ 
             ARU +
             s(Anthropophony, by = ARU) + 
             s(Biophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = scat(link = "identity"),
           method = 'REML')
AIC(mod5, mod6)
summary(mod6)
gam.check(mod6) # better diagnostics

# doesnt work with slope offsets (June 13)
#pred_plot(data_df = temp_df, model_fit = mod6, index_name = "logR")

###########################################
################# rugo ####################
i = "rugo"
temp_y = data.frame(site = indices_df$site, rugo = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU))

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
gam.check(mod2) # much better error structure in QQ and variance in residuals.

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

###########################################
################# sfm #####################
i = "sfm"
temp_y = data.frame(site = indices_df$site, sfm = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU))

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

# visualize pairs of corr
# temp_df %>%
#   GGally::ggpairs(aes(alpha = 0.05))

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
AIC(mod2, mod3)

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
AIC(mod2, mod4) #deltaAIC >3 and mod4 is more parsimonious

pred_plot(data_df = temp_df, model_fit = mod4, index_name = i)

###########################################
################# zcr_mean ################
i = "zcr_mean"
temp_y = data.frame(site = indices_df$site, zcr_mean = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU))

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

# visualize pairs of corr
# temp_df %>%
#   GGally::ggpairs(aes(alpha = 0.05))

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

# Quiet has complex fit/higher basis needed
mod3 = gam(zcr_mean ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet, k = 25) + 
             s(Interference) +
             s(ARU, bs = "re", k = 2),
           data = temp_df,
           family = gaussian,
           method = 'REML')
gam.check(mod3) # QQ plot deviates a bit more than without RE but model is highly recommended over no RE
summary(mod3)
AIC(mod2, mod3) 


draw(mod2, scales = "fixed")
pred_plot(data_df = temp_df, model_fit = mod2, index_name = i)


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

