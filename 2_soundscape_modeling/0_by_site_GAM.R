# GAM analysis of acoustic indices

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
source('/paper-AcousticIndices/utility_fxs.R')

###########################################
####### read in data and data prep ########
wd = '/paper-AcousticIndices/results/'

# Read in ABGQI classified site average data and organize
abgqi_df = fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI.csv')) %>%
  mutate(ARU = substr(site, 4,5),
         ARUdevice = substr(site, 4,8)) %>%
  select(-contains('var')) %>%
  select(site, wavs, contains('mean'), ARU, -Unidentified_mean) %>%
  rename(Anthropophony = Anthropophony_mean, Biophony = Biophony_mean, Geophony = Geophony_mean,
         Quiet = Quiet_mean, Interference = Interference_mean)
  
# Site average acoustic index csv
indices_df = fread(paste0(wd, 'acoustic_indices_aggregation/averages/site_avg_acoustic_indices.csv'))

# Annual summary
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
  filter(wavs >= (144 * 2.5)) # rm 44 sites

# sites with static
# s2lam050_190411, s2lam049_210501 : a lot of static but still has other signal
error_sites = read.csv('/paper1-AcousticIndices/data/identified_problem_sites_witherrors.csv')
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

# scale covariates (ABGQI with 0-1 min max)
mod_df = mod_df %>%
  mutate(Anthropophony = min_max_norm(Anthropophony),
         Biophony      = min_max_norm(Biophony),
         Geophony      = min_max_norm(Geophony),
         Quiet         = min_max_norm(Quiet),
         Interference  = min_max_norm(Interference),
         ARU = factor(ARU),
         logwavs = log(wavs)) 

# save list of sites used in analyses
# sites = temp_df$site
# write.csv(data.frame(sites), 
#           paste0(wd, 'sites_used_in_GAMs_21July22.csv'), 
#           row.names = F)

# list to store model slopes that are saved
model_slopes = list()
final_models = list()



# Note on GAMs:
# - every acoustic index is manually modeled using backward selection using ARU, n Recordings, ABGQI (see Supplementary Material for details) 
# - Final GAM fits are stored in a model list and saved as an RData object at the end of script
# - Index GAMs are fit in alphabetical order below

###########################################
################# ACI #####################
i = "ACI"
temp_y = data.frame(site = indices_df$site, ACI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df)

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
# bring ACI down to just above 0 (for log transform)
temp_df$ACI_0 = temp_df$ACI - min(temp_df$ACI) + 1e-7

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

# remove extreme residual
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
remod1$aic
mod2$aic

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
mod2$aic
mod3$aic
gam.check(mod3, rep = 500)
draw(mod3, scales = "fixed") 
appraise(mod3, method = "simulate")

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod3

# visualize predictions vs observed
pred_plot(data_df = temp_df, model_fit = mod3, index_name = "ACI_0")

# clean up models
rm(temp_df); rm(i); rm(mod1); rm(mod2); rm(mod3); rm(remod1)

###########################################
################# ADI #####################
i = "ADI"
temp_y = data.frame(site = indices_df$site, ADI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df)

# visualize data
temp_df %>%
  gather(soundType, value, -ADI, -ARU, -wavs, -site) %>%
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
gam.check(remod2, rep = 500) # Okay. Some odd fit issues in resids plot but are still roughly normal. Two low QQ tail deviations


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
mod2$aic
mod3$aic
gam.check(mod3, rep = 500)

draw(mod3, scales = "fixed")

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod3

# visualize predictions vs observed
pred_plot(data_df = temp_df, model_fit = mod3, index_name = "ADInorm")

# clean up models
rm(temp_df); rm(i); rm(mod1); rm(mod2); rm(mod3); rm(remod2)

###########################################
################# AEI #####################
i = "AEI"
temp_y = data.frame(site = indices_df$site, AEI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df)

# visualize data
temp_df %>%
  gather(soundType, value, -AEI, -ARU, -wavs, -site) %>%
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
mod1$aic
mod2$aic
gam.check(mod2)

draw(mod2, scales = "fixed")

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod2

# clean up models
rm(temp_df); rm(i); rm(mod1); rm(mod2)

###########################################
################# BI ######################
i = "BI"
temp_y = data.frame(site = indices_df$site, BI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) 

# visualize data
temp_df %>%
  gather(soundType, value, -BI, -ARU, -wavs, -site) %>%
  ggplot(aes(x = value, y = BI, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")

# both ARUs have similar distribution - start with gamma
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
gam.check(mod2)
mod1$aic
mod2$aic

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
mod2$aic
mod3$aic
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
rm(temp_df); rm(i); rm(mod1); rm(mod2); rm(mod3)

###########################################
################# H #######################
i = "H"
temp_y = data.frame(site = indices_df$site, H = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) 

# visualize data
temp_df %>%
  gather(soundType, value, -H, -ARU, -wavs, -site) %>%
  ggplot(aes(x = value, y = H, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = H, fill = ARU)) +
  geom_histogram()
summary(temp_df$H)


# fit a model bound to 0-1
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
mod1$aic
mod2$aic


# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod2

pred_plot(data_df = temp_df, model_fit = mod2, index_name = i)

# clean up models
rm(temp_df); rm(i); rm(mod1); rm(mod2)

###########################################
################# Hs ######################
i = "Hs"
temp_y = data.frame(site = indices_df$site, Hs = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) 

# visualize data
temp_df %>%
  gather(soundType, value, -Hs, -ARU, -wavs, -site) %>%
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
mod1$aic
mod2$aic

pred_plot(data_df = temp_df, model_fit = mod2, index_name = i)

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod2

# clean up models
rm(temp_df); rm(i); rm(mod1); rm(mod2)

###########################################
################# Ht ######################
i = "Ht"
temp_y = data.frame(site = indices_df$site, Ht = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) 

# visualize data
temp_df %>%
  gather(soundType, value, -Ht, -ARU, -wavs, -site) %>%
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
gam.check(mod1) # good diags - maybe pronounced residual shoulders

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
mod1$aic
mod2$aic


# Model 2
pred_plot(data_df = temp_df, model_fit = mod2, index_name = i)

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod2

# Clean up
rm(temp_df); rm(i); rm(mod1); rm(mod2)

###########################################
################# M ######################
i = "M"
temp_y = data.frame(site = indices_df$site, M = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df)

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -M, -ARU, -wavs, -site) %>%
  ggplot(aes(x = value, y = M, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")
ggplot(temp_df, aes(x = M, fill = ARU)) +
  geom_histogram()
summary(temp_df$M)

# Extreme high values are in M
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
mod3$aic
mod4$aic

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
mod4$aic
mod5$aic # suggests keeping Quiet

pred_plot(data_df = temp_df_remod, model_fit = mod4, index_name = "logM")

# summarized slopes
ci = slope_summary(mod4)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod4

# Clean up
rm(temp_df); rm(i); rm(mod1); rm(mod2); rm(mod3); rm(mod4); rm(mod5); rm(remod2);


###########################################
################# NDSI ####################
i = "NDSI"
temp_y = data.frame(site = indices_df$site, NDSI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df)

# visualize data
temp_df %>%
  gather(soundType, value, -NDSI, -ARU, -wavs, -site) %>%
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
gam.check(mod1) 


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
mod1$aic
mod2$aic

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
mod2$aic
mod3$aic

draw(mod3, scales = "fixed", )
vis.gam(mod3, theta = 65)
gam.vcomp(mod3)

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod3

# visualize predictions vs observed
pred_plot(data_df = temp_df, model_fit = mod3, index_name = "beta_NDSI")

# Clean up
rm(temp_df); rm(i); rm(mod1); rm(mod2); rm(mod3)

###########################################
################# NDSI_A ##################
i = "NDSI_A"
temp_y = data.frame(site = indices_df$site, NDSI_A = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df)

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -NDSI_A, -ARU, -wavs, -site) %>%
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
mod1$aic
mod2$aic

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
mod2$aic
mod3$aic

pred_plot(data_df = temp_df, model_fit = mod3, index_name = i)

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod3

# clean up 
rm(temp_df); rm(i); rm(mod1); rm(mod2); rm(mod3); 

###########################################
################# NDSI_B ##################
i = "NDSI_B"
temp_y = data.frame(site = indices_df$site, NDSI_B = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) 

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -NDSI_B, -ARU, -wavs, -site) %>%
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
remod2$aic
mod3$aic

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
mod3$aic
mod4$aic

pred_plot(data_df = temp_df, model_fit = mod4, index_name = "normNDSI_B")

# summarized slopes
ci = slope_summary(mod4)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod4

# Clean up
rm(temp_df); rm(i); rm(mod1); rm(mod2); rm(mod3); rm(mod4); rm(remod2)

###########################################
################# R #######################
i = "R"
temp_y = data.frame(site = indices_df$site, R = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) 

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -R, -ARU, -wavs, -site) %>%
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
summary(mod3)
gam.check(mod3)
mod2$aic
mod3$aic

# drop logwavs (p = 0.15928)
mod4 = gam(logR ~ 
             ARU +
             s(Anthropophony, k = 5) + 
             s(Biophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = temp_df,
           method = 'ML')
summary(mod4)
gam.check(mod4)
mod3$aic
mod4$aic

pred_plot(data_df = temp_df, model_fit = mod4, index_name = 'logR')

# summarized slopes
ci = slope_summary(mod4)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod4

# Clean up
rm(temp_df); rm(i); rm(mod1); rm(mod2); rm(mod3); rm(mod4)

###########################################
################# rugo ####################
i = "rugo"
temp_y = data.frame(site = indices_df$site, rugo = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) 

# visualize data
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -rugo, -ARU, -wavs, -site) %>%
  ggplot(aes(x = value, y = rugo, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = 'free')
ggplot(temp_df, aes(x = rugo, fill = ARU)) +
  geom_histogram()
summary(temp_df$rugo) # positive continuous


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
gam.check(mod2)
remod1$aic
mod2$aic

pred_plot(data_df = temp_df, model_fit = mod2, index_name = i)

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod2

# Clean up
rm(temp_df); rm(i); rm(mod1); rm(mod2); rm(remod1)

###########################################
################# sfm #####################
i = "sfm"
temp_y = data.frame(site = indices_df$site, sfm = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df)

# visualize data - definitely different y-intercepts for ARU
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -sfm, -ARU, -wavs, -site) %>%
  ggplot(aes(x = value, y = sfm, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")

# approximately normal 0-1 bound or positive continuous
ggplot(temp_df, aes(x = sfm, fill = ARU)) +
  geom_histogram()
summary(temp_df$sfm)


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
remod1$aic
mod2$aic

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
mod2$aic
mod3$aic # recommended to keep geophony

pred_plot(data_df = temp_df, model_fit = mod2, index_name = i)

# summarized slopes
ci = slope_summary(mod2)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod2

# Clean up
rm(temp_df); rm(i); rm(mod1); rm(mod2); rm(mod3); rm(remod1)

###########################################
################# zcr_mean ################
i = "zcr_mean"
temp_y = data.frame(site = indices_df$site, zcr_mean = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) 

# visualize data - definitely different y-intercepts for ARU
temp_df %>%
  select(-site, -wavs) %>%
  gather(soundType, value, -zcr_mean, -ARU, -wavs, -site) %>%
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
gam.check(mod2) # hard to tell if diags are better - one extreme value

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
remod2$aic
mod3$aic

draw(mod3, scales = "fixed")
pred_plot(data_df = temp_df, model_fit = mod3, index_name = "normZCR")

# summarized slopes
ci = slope_summary(mod3)
model_slopes[[i]] = cbind(index = rep(i, times = nrow(ci)), ci)
final_models[[i]] = mod3

# Clean up
rm(temp_df); rm(i); rm(mod1); rm(mod2); rm(mod3); rm(remod2);

###########################################
######### SAVE OBJECTS ####################
saveRDS(final_models, file = '/models/gam_model_objects.RData')
saveRDS(model_slopes, file = '/models/gam_model_slopes.RData')

