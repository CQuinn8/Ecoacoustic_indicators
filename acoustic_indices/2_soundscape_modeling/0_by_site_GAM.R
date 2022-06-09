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
# ACI: (June9) complete: a few small deviations in qq plot otherwise okay
# ADI: (June9) complete: used beta distrubition on min max normalized values [0-1]
# AEI: (June9) complete: relatively simple but sufficient model
# BI: (June9) complete - possible model visual started too
# H: (June9) complete
# Hs: (June9) complete: model may be too simple -revisit
# Ht:
# M:
# NDSI_A:
# NDSI_B:
# R:
# rugo:
# sfm:
# zcr_mean:


# TODO:
# finsish initial model runs
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

# Clean sites with anamolous values
# 1. Top 0.1% of Sites with significant interference (e.g., Int > 4sd(Int))
# 2. Sites with any value == 1
# 3. Zero averages -> median value
# Interference
int_sd = sd(abgqi_df$Interference)
mod_df = abgqi_df %>%
  filter(Interference <= 4 * int_sd) # 20 sites

# only one site with any value == 1 (quiet)
mod_df = mod_df %>%
  filter(Quiet < 1)

# account for zero values
mod_df[mod_df == 0] = NA

# Which rows contain zeroes
zeros_df = mod_df[!complete.cases(mod_df),] # n = 24: 22 are quiet, 1 Anthro & Geo, 1 Int

# calculate medians
medians = mod_df %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)))

# input medians for 0 (now NA) values
mod_df = mod_df %>% 
  mutate(across(where(is.numeric), .fns = ~ifelse(is.na(.x), median(.x, na.rm=TRUE), .x)))

# Number of minutes in model data
sum(mod_df$wavs) #731634

###########################################
################# NDSI ####################
i = "NDSI"
temp_y = data.frame(site = indices_df$site, NDSI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU)) %>%
   select(-c(site))

# visualize data
temp_df %>%
  gather(soundType, value, -NDSI, -ARU) %>%
  ggplot(aes(x = value, y = NDSI, colour = ARU)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "gam") +
    facet_wrap(~soundType)
  
# add transformed NDSI
temp_df$posNDSI = temp_df$NDSI + 1
temp_df$beta_NDSI = (temp_df$NDSI + 1)/2
temp_df$sq_beta_NDSI = (temp_df$beta_NDSI)^2
hist(temp_df$beta_NDSI)

# visualize pairs of corr
temp_df %>%
  GGally::ggpairs(aes(alpha = 0.05), progress = FALSE)

###### TESTING using: https://www.youtube.com/watch?v=Ukfvd8akfco
mod1 = gam(NDSI ~ 
              Anthropophony + 
              Biophony + 
              Geophony + 
              Quiet + 
              Interference,
            data = temp_df)
summary(mod1)

# model using beta transformed NDSI
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

AIC(mod1, mod2)

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
variance_comp(mod3)


# Add random slopes
mod4 = gam(beta_NDSI ~ 
             ARU + # must include base factor to interpret offset of "by"
             s(Anthropophony, by = ARU, bs = 're') + 
             s(Biophony, by = ARU, bs = 're') + 
             s(Geophony, by = ARU, bs = 're') + 
             s(Quiet, by = ARU, bs = 're') + 
             s(Interference, by = ARU, bs = 're') +
             s(ARU, bs = 're'),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod4)

AIC(mod3, mod4)

# Model 3 appears more parsimoniuous
# remove geophony
mod5 = gam(beta_NDSI ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're'), # random intercept
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod5)
AIC(mod3, mod5)


# model 3 is best option
mod = mod3
draw(mod) # nice visualization tool from gratia for partial effects

k.check(mod) 

#appraise(mod2, method = "simulate") 
# QQ should follow 1:1 wihtin the reference band (do data follow assumptions of model well) some do not
# residuals should be approx gaussian 
par(mfrow = c(2, 2))
vis.gam(mod, theta = 65)

gam.check(mod, rep = 500) # for rep = 500 does the QQ roughly fall within confidence interval? Not at the upper tail
                           # histogram is roughly gaussian, Resp vs fitted is roughly 1:1, resids vs linear pred may have decreasing variance

mod6 = gam(beta_NDSI ~ 
             s(Anthropophony, k =100) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet, k =100) + 
             s(Interference) +
             s(ARU, bs = 're'), # random intercept
           data = temp_df,
           family = betar(),
           method = 'REML')

summary(mod6) # model summary
AIC(mod3, mod6)

draw(mod6)

gam.check(mod6)


# credible intervals 
pred = predict.gam(mod, type = "lpmatrix")


###########################################
################# ACI #####################
i = "ACI"
temp_y = data.frame(site = indices_df$site, ACI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  mutate(ARU = factor(ARU)) %>%
  select(-c(site,wavs))

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

# remove one extremely high value
outlier(temp_df$ACI)
temp_df = temp_df %>%
  filter(ACI < outlier(ACI))
# faraway::halfnorm(temp_df$ACI)

# visualize pairs of corr
temp_df %>%
  GGally::ggpairs(aes(alpha = 0.05))

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
draw(mod3, scales = "fixed") 
appraise(mod3, method = "simulate") # only works with gaussian data
gam.check(mod3, rep = 500)

# increase anthropophony basis functions per gam.check
mod4 = gam(ACI_0 ~ 
             s(Anthropophony, k = 25) + 
             s(Biophony) + 
             s(Geophony) + 
             s(Quiet) + 
             s(Interference) +
             s(ARU, bs = 're', k = 2),
           data = temp_df,
           family = Gamma(link = "log"),
           method = 'REML')
summary(mod4)
gam.check(mod4) # Anthro is >0.05 and edf is << than k'
AIC(mod3, mod4)
draw(mod4, scales = "fixed")


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

mean(temp_df$ADI);var(temp_df$ADI)

# add transformed 
temp_df$ADI3 = temp_df$ADI^3
temp_df$ADInorm = min_max_norm(temp_df$ADI, min_x = min(temp_df$ADI), max_x = max(temp_df$ADI)) # bring to 0-1 scale

# remove one extremely low value
outlier(temp_df$ADI)
temp_df = temp_df %>%
  filter(ADI > outlier(ADI))

# visualize pairs of corr
temp_df %>%
  GGally::ggpairs(aes(alpha = 0.05))


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
draw(mod2) 
appraise(mod2, method = "simulate") # better behavior but still heterosk in resids
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
draw(mod3, scales = "fixed") 
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
draw(mod4, scales = "fixed") # random effect reduces the complexity of splines somewhat

# check for response scale compared to raw data 
mean(fitted(mod1));mean(temp_df$ADI) # approximates exactly
mean(fitted(mod3)); mean(fitted(mod4)); mean(temp_df$ADInorm) # is within a reasonable margin ( - 0.08)


###########################################
################# ADI #####################
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
ggplot(temp_df, aes(x = sqrtAEI, fill = ARU)) +
  geom_histogram()

# add transformed 
temp_df$sqrtAEI = sqrt(temp_df$AEI) # approximates normal distribution


# visualize pairs of corr
temp_df %>%
  GGally::ggpairs(aes(alpha = 0.05))


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
par(mfrow = c(2, 2))
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

# visualize pairs of corr
temp_df %>%
  GGally::ggpairs(aes(alpha = 0.05))

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

# PRED VISUAL
# make the prediction, add this and a column of standard errors to the prediction data.frame.
# COULD USE THIS VISUAL AND SPECIFY EACH K VALUE ACCORDING TO MODEL
pred = data.frame(
  Anthropophony = temp_df$Anthropophony,
  Biophony = temp_df$Biophony,
  Geophony = temp_df$Geophony,  
  Quiet = temp_df$Quiet,
  Interference = temp_df$Interference,
  BI = temp_df$BI,
  predict(mod4, type = "response", se = TRUE)) %>%
  gather(variable, value, -fit, -se.fit, -BI)

(ggAnthro = pred %>%
  filter(variable == "Biophony") %>%
  ggplot(aes(x = value)) +
  geom_point(aes(y = BI), size = 1, alpha = 0.5) +
  geom_line(aes(y = fit), alpha = 0.3) + # fit values
  geom_smooth(aes(y = fit), colour = "red", alpha = 0.3, method = "gam", formula = y ~ s(x, k = 5))) # smoothed fit

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

# remove one extremely low value
outlier(temp_df$H)
temp_df = temp_df %>%
  filter(ADI > outlier(ADI))

# visualize pairs of corr
temp_df %>%
  GGally::ggpairs(aes(alpha = 0.05))


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
AIC(mod1, mod2)
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

# visualize pairs of corr
temp_df %>%
  GGally::ggpairs(aes(alpha = 0.05))


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
AIC(mod1, mod2)
gam.check(mod2) # better diagnostics, still some influential points though in QQ lower end, low shoulders in residuals histogram

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
gam.check(mod3) # slightly worse QQ patterns, residuals have possible taper pattern, histogram still has sharp peak
draw(mod3, scales = "fixed")


####

















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
# # set up pred df
# temp_df$y_pred = predict(fit, newdata = temp_df)
# 
# temp_df %>%
#   gather(covariate, value, -NDSI, -y_pred) %>%
#   ggplot(aes(x = value, y = y_pred)) +
#   geom_point() +
#   geom_point(aes(y = NDSI), colour = 'red', alpha = 0.2) +
#   facet_wrap(~covariate, scales = 'free_x')
# 
# temp_df %>%
#   mutate(y_diff = NDSI - y_pred) %>%
#   gather(covariate, value, -NDSI, -y_pred, -y_diff) %>%
#   ggplot(aes(x = value, y = y_diff)) +
#   geom_point() +
#   facet_wrap(~covariate, scales = 'free_x')
# 
# ###########################################
# # ACI
# i = responses[1]
# temp_y = data.frame(site = indices_df$site, ACI = indices_df[[i]])
# temp_df = merge(x = temp_y, y = mod_df) %>% 
#   select(-c(site,wavs))
# 
# ggplot(temp_df, aes(x = 1/(ACI)^2)) + geom_histogram()
# ggplot(temp_df, aes(x = log(ACI))) + geom_histogram()
# 
# # visualize pairs of corr
# df %>%
#   GGally::ggpairs(aes(alpha = 0.05))
# 
# temp_df %>% 
#   mutate(y = log(y)) %>%
#   response_exam(i)
# 
# # FIT GAM
# fit1 = gam(ACI ~ s(Anthropophony_mean) + s(Anthropophony_mean) + s(Biophony_mean) + s(Geophony_mean) + 
#              s(Quiet_mean) + s(Interference_mean) + s(Unidentified_mean), 
#            data = temp_df)
# summary(fit1)
# 
# plot.gam(fit1, residuals = TRUE, # residuals
#          shade = TRUE, # 95% CI for the mean shape of effects
#          shade.col = "lightblue", seWithMean = TRUE, pages = 1, pch =1, cex = 1, all.terms = TRUE)
# 
# plot(fit1, pages = 1, all.terms = TRUE) # partial effects plots
# 
# # small p-values indicate residuals are not randomly distributed (i.e., not enough basis fxs)
# gam.check(fit1)
# # if p is significant, try increasing k: gam(y ~ s(x1, k = 12))
# 
# # FIT 2
# # geophony is weakly significant, increase k
# fit2 = gam(log(ACI) ~ s(Anthropophony_mean) + s(Anthropophony_mean) + s(Biophony_mean) + s(Geophony_mean, k = 12) + 
#              s(Quiet_mean) + s(Interference_mean) + s(Unidentified_mean), 
#            data = temp_df)
# summary(fit2)
# gam.check(fit2)
# 
# # values over ~0.8 for worst may mean issues with collinearity
# concurvity(fit2, full = TRUE)
# concurvity(fit2, full = FALSE) # use to id where the relationships are UNID and Int is highest
# 
# y_pred = predict(fit2, temp_df[,2:7])
# 
# 
# # FIT 3 TRYING WITH TRANSFORMED ACI
# fit3 = gam(log(ACI) ~ s(log(Anthropophony_mean)) + s(log(Anthropophony_mean)) + s(log(Biophony_mean)) + 
#              s(log(Geophony_mean)) + s(log(Quiet_mean)) + s(log(Interference_mean)) + s(log(Unidentified_mean)), 
#            data = temp_df)
# summary(fit3)
# gam.check(fit3)
# plot(fit3, pages = 1, all.terms = TRUE) # partial effects plots
# autoplot(fit3)
