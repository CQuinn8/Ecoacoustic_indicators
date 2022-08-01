# Use by-site CNN counts of pres/abs
# 1. Derive Spp Richness (e.g., cnt >= 3 == presence)
# 2. GAM: SppRichness ~ ABGQI
# 3. AcInd: SppRichness ~ [indices]

# data libs
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)


# modeling libs
library(mgcv)
library(gratia)
source('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/code/paper1-AcousticIndices/utility_fxs.R')

wd = '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/'

# bird species pres/abs by site
bird_df = fread('G:/Shared drives/NASA  S2L/Paper Development/Sonoma SDM update/bird_cnn_predictions/sonoma_s2l_predictions_ROImaxF05_220707_summarized.csv')

# reduce sites used in ABGQI analyses
sites = fread("//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/sites_used_in_GAMs_21July22.csv")

# subsetting robust ABGQI sites results in 1221 -> 1182 sites
bird_df = bird_df[bird_df$site %in% sites$sites, ]

# count number of occ per spp
spp_classifications = bird_df %>%
  select(-Richness, -site) %>%
  colSums()

# default n > 3 observations for species to be present
site_spp_rich = bird_df %>%
  select(-Richness) %>%
  gather(spp, count, -site) %>% # convert df into long format
  mutate(presAbs = ifelse(count >= 3, 1, 0)) %>% # mutate new column into pres/abs if count > 3
  group_by(site) %>% # group by site 
  summarise(site_richness = sum(presAbs)) # count the number of species that are present

# visualize distribution of spp richness
site_spp_rich %>%
  ggplot(aes(x = site_richness)) +
  geom_histogram()


# read in ABGQI and Acoustic Indices
abgqi_df = fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI.csv')) %>%
  mutate(ARU = substr(site, 4,5),
         ARUdevice = substr(site, 4,8)) %>%
  select(-contains('var')) %>%
  select(site, wavs, contains('mean'), ARU, -Unidentified_mean) %>%
  rename(Anthropophony = Anthropophony_mean, Biophony = Biophony_mean, Geophony = Geophony_mean,
         Quiet = Quiet_mean, Interference = Interference_mean)
indices_df = fread(paste0(wd, 'acoustic_indices_aggregation/averages/site_avg_acoustic_indices.csv'))


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

# scale covariates ABGQI
mod_df = mod_df %>%
  mutate(Anthropophony = min_max_norm(Anthropophony),
         Biophony      = min_max_norm(Biophony),
         Geophony      = min_max_norm(Geophony),
         Quiet         = min_max_norm(Quiet),
         Interference  = min_max_norm(Interference))
mod_df = merge(site_spp_rich, mod_df, by = 'site')

mean(mod_df$site_richness)
sd(mod_df$site_richness)

mod_df = mod_df %>%
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) %>%
  select(-c(site))

#####################################
# SPP RICH ~ SOUNDSCAPE SOURCES
# visualize data
mod_df %>%
  gather(soundType, value, -site_richness, -ARU) %>%
  ggplot(aes(x = value, y = site_richness, colour = ARU)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "gam") +
  facet_wrap(~soundType, scales = "free")

# start with full model with poisson
mod1 = gam(site_richness ~ 
             ARU +
             logwavs +
             s(Anthropophony, k = 5) + 
             s(Geophony, k = 5) + 
             s(Quiet, k = 5) + 
             s(Interference, k = 5),
           data = mod_df,
           family = poisson(),
           method = 'ML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # looks good
draw(mod1, scales = 'fixed')

# summarized slopes
ci = slope_summary(mod1)
model_slopes = cbind(index = rep('Bird Spp. Richness', times = nrow(ci)), ci)

# visualize predictions vs observed
pred_plot(data_df = mod_df, model_fit = mod1, index_name = "site_richness")

# Visualize slope objects 
# calculate all slope values
slope_df = model_slopes

# use middle 99% of x values to plot slopes to avoid extreme tail behavior
abgqi_limits = mod_df %>%
  select(Anthropophony, Biophony, Geophony, Quiet, Interference) %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(lower = max(0, mean(value) - 3 * sd(value)), # set to zero or sd value if >0
            upper = min(1, mean(value) + 3 * sd(value))) # set to 1 or sd value if < 1

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
  group_by(var) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)))

temp_min = floor(min(slope_df_filtered$lower))
temp_max = ceiling(max(slope_df_filtered$upper))
(gg = slope_df_filtered %>%
    group_by(var) %>%
    mutate(var = factor(var, levels = c('Geophony', 'Anthropophony', 'Interference', 'Quiet'))) %>%
    ggplot(aes(x = var, y = derivative)) +
    geom_boxplot(width = 0.25, outlier.shape = NA) +
    #geom_jitter(width = 0.05, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ylim(c(temp_min, temp_max)) + # min max in slope df
    ylab("Slope") +
    xlab("Soundscape Source") +
    coord_flip() +
    theme_bw())
annotate_figure(allplots, top = text_grob("95% CI for first derivative (slopes) of GAM partial effects.
Note: should be interpreted alongside partial effects plots.", size = 14, face = "bold"),
                bottom = text_grob("Values reflect middle 99% of covariate range based on erroneous tail behavior.",
                                   size = 10)) 

#############################################
# BIOPHONY ~ SPP RICH 
# visualize data
mod_df %>%
  select(Biophony, site_richness, ARU) %>%
  ggplot(aes(y = Biophony, x = site_richness, colour = ARU)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "gam")
hist(mod_df$Biophony)
hist(sqrt(mod_df$Biophony))

mod_df$sqrtBiophony = sqrt(mod_df$Biophony)
mod_df$cuberootBiophony = (mod_df$Biophony)^(1/3)
mod_df$logBiophony = log(mod_df$Biophony+1e-7)
hist(mod_df$logBiophony) # one extreme outlier that will impact model fit [886,]
mod_df = mod_df[-886,]

modBio1 = gam(cuberootBiophony ~ 
                ARU +
                logwavs +
                s(site_richness, k = 5),
             data = mod_df,
             family = betar(),
             method = 'ML')
summary(modBio1)
par(mfrow = c(2,2))
gam.check(modBio1) # some tail behavior in QQ and a right skewed histogram - check on this value

# check on the extreme residuals
mod_df$resids = data.frame(residuals.gam(modBio1, type = 'deviance'))
mod_df_remod = mod_df %>% # [-639,] 
  filter(resids < 8) # high Biophony (>95%) with low spp richness

remodBio1 = gam(cuberootBiophony ~ 
                ARU +
                logwavs +
                s(site_richness, k = 5),
              data = mod_df_remod,
              family = betar(),
              method = 'ML')
summary(remodBio1)
gam.check(remodBio1) # qq plot still showing right skewed data but histogram is better

draw(remodBio1)

# RMSE
sqrt(mean((mod_df_remod$Biophony - remodBio1$fitted.values^3)^2))

# Visualize model predictions and actual
# predict based on original values
mod_df_remod$Biophony_pred = predict(remodBio1, newdata = mod_df_remod, type = "response")^3

mod_df_remod %>%
  select(site_richness, Biophony, Biophony_pred, ARU) %>%
  gather(covariate, value, -Biophony, -Biophony_pred, -ARU) %>%
  ggplot(aes(x = value, y = toLogit(Biophony_pred))) +
    geom_point(aes(x = value, y = toLogit(Biophony), colour = ARU), alpha = 0.5) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = 'gam', colour = 'black', alpha = 0.4) +
    ggtitle("Biophony: black = predicted values") +
    labs(y = 'Percent Biophony', x = 'Bird Species Richness') +
    theme_bw()

pred_plot(data_df = mod_df_remod, model_fit = remodBio1, index_name = "cuberootBiophony")

#############################################
# SPP RICH ~ ACOUSTIC INDEX ANALYSIS
mod_df_indices = merge(indices_df, site_spp_rich, by = 'site') 
mod_df_indices = merge(mod_df_indices, abgqi_df_no_error, by = 'site') %>%
  mutate(logwavs = log(wavs),
         ARU = factor(ARU)) %>%
  select(-YYYY, -site, -wavs, -Anthropophony, -Biophony, -Geophony, -Quiet, -Interference)
  

mod1 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+ 
             s(Hs, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(sfm, k = 5)+ 
             s(rugo, k = 5)+ 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod1)
gam.check(mod1, type = 'response') # correct family and error structure

# drop rugo (p = 0.9233)
mod2 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+ 
             s(Hs, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(sfm, k = 5)+ 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod2)

# drop zcr (p = 0.4670)
mod3 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+ 
             s(Hs, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(sfm, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod3)


# drop ARU (p = 0.33)
mod4 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+ 
             s(Hs, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(sfm, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod4)

# drop Hs (p = 0.2042)
mod5 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+  
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(sfm, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod5)

# try dropping R (p = 0.051070)
mod6 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+  
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(sfm, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')

AIC(mod5, mod6) # AIC suggests keeping R in model

# mod
summary(mod5)
par(mfrow = c(2,2))
gam.check(mod5, type = 'response')
draw(mod5, scales = 'fixed')


# summarized slopes
ci = slope_summary(mod5)
model_slopes = cbind(index = rep('Bird Spp. Richness', times = nrow(ci)), ci)

# Visualize slope objects 
# calculate all slope values
slope_df = model_slopes

# Use middle 95% of data for each index domain to eliminate spurious tail behavior in plots
slope_df_filtered = slope_df %>%
  select(-smooth) %>%
  group_by(var) %>%
  filter(data >= quantile(data, 0.025) & data <= quantile(data, 0.975))

temp_min = floor(min(slope_df_filtered$lower))
temp_max = ceiling(max(slope_df_filtered$upper))
(gg = slope_df_filtered %>%
    group_by(var) %>%
    ggplot(aes(x = var, y = derivative)) +
      geom_boxplot(width = 0.25, outlier.shape = NA) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ylim(c(temp_min, temp_max)) + # min max in slope df
      ylab("Slope") +
      xlab("Acoustic Index") +
      coord_flip() +
      theme_bw())
annotate_figure(gg, top = text_grob("95% CI for first derivative (slopes) of GAM partial effects.
Note: should be interpreted alongside partial effects plots.", size = 14, face = "bold"),
                bottom = text_grob("Values reflect middle 99% of covariate range based on erroneous tail behavior.",
                                   size = 10)) 

