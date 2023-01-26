
# data libs
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)

# modeling libs
library(mgcv)
library(gratia)
library(lmtest)

source('utility_fxs.R')

###################################
# DATA IMPORT
# bird species pres/abs by site
bird_df = fread('data/sonoma_s2l_predictions_SoundscapeMaxF05_230118_summarized-dawn_4am-12pm.csv')

# reduce sites used in ABGQI analyses
sites = fread("data/sites_used_in_GAMs.csv")

# ABGQI
abgqi_df = fread('data/site_ABGQI_dawn_4am-12pm.csv') %>%
  mutate(ARU = substr(site, 4,5)) %>%
  select(-Unidentified)

# Acoustic Indices
indices_df = fread('data/site_acoustic_indices_dawn_4am-12pm.csv')

# sites with static
error_sites = read.csv('data/identified_problem_sites_witherrors.csv')
error_sites = paste0(error_sites$SiteID, collapse = '|')

###################################
# Prep data and summarize
# subsetting robust ABGQI sites results in 1226 -> 1185 sites
bird_df = bird_df[bird_df$site %in% sites$sites, ]
abgqi_df_no_error = abgqi_df[abgqi_df$site %in% sites$sites,]
indices_df = indices_df[indices_df$site %in% sites$sites, ]

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

# Merge dataframes dropping the 10 sites with no bird spp observations
mod_df = abgqi_df_no_error %>%
  left_join(y = site_spp_rich, by = 'site') %>%
  drop_na()

# Stats on species richness and biophony
mean(mod_df$site_richness)
sd(mod_df$site_richness)
mean(mod_df$Biophony) # 0.3501
sd(mod_df$Biophony) # 0.1918

mod_df = mod_df %>%
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) 

mod_df$sqrtBiophony = sqrt(mod_df$Biophony)
mod_df$cuberootBiophony = (mod_df$Biophony)^(1/3)
mod_df$logBiophony = log(mod_df$Biophony+1e-7)

#############################################
# Eq 2 Dawn: BIOPHONY ~ SPP RICH 
# visualize data
mod_df %>%
  select(Biophony, site_richness, ARU) %>%
  ggplot(aes(y = Biophony, x = site_richness, colour = ARU)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "gam")

hist(mod_df$Biophony)
hist(mod_df$sqrtBiophony)
hist(mod_df$cuberootBiophony)
hist(mod_df$logBiophony)

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
concurvity(modBio1)

# RMSE
sqrt(mean((mod_df$Biophony - modBio1$fitted.values^3)^2))

termplot(modBio1)

# Visualize model predictions and actual
# predict based on original values
mod_df$Biophony_pred = predict(modBio1, newdata = mod_df, type = "response")^3

# Grouped by ARU
(gg <- mod_df %>%
  select(site_richness, Biophony, Biophony_pred, ARU) %>%
  gather(covariate, value, -Biophony, -Biophony_pred, -ARU) %>%
  ggplot(aes(x = value, y = Biophony_pred * 100)) +
    geom_point(aes(x = value, y = Biophony * 100, colour = ARU), alpha = 0.5) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = 'gam', colour = 'black', alpha = 0.4) +
    # ggtitle("Biophony: black = predicted values") +
    labs(y = 'Percent Biophony', x = 'Bird Species Richness') +
    scale_color_hue(labels = c("Audiomoth", "LG"), ) +
    theme_bw() +
    guides(colour = guide_legend(title = 'Observed values')) +
    theme(legend.position = "top"))


#############################################
# Eq 3 Dawn: SPP RICH ~ Biophony
mod_df_bio = mod_df %>%
  select(Biophony, site_richness, logwavs, ARU) %>%
  mutate(Biophony = min_max_norm(Biophony))

mod1 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5),
           data = mod_df_bio,
           family = poisson(),
           method = 'ML')

summary(mod1)
gam.check(mod1)
draw(mod1)
concurvity(mod1)

# RMSE
sqrt(mean((mod_df_bio$site_richness - mod1$fitted.values)^2))

# rRMSE 
sqrt(mean((mod_df_bio$site_richness - mod1$fitted.values)^2)) / 
  diff(range(mod_df_bio$site_richness))


mod_df_bio$spp_rich_pred = predict(mod1, newdata = mod_df_bio, type = "response")
pred_plot(mod_df_bio, mod1, "site_richness")


#############################################
# Eq 4 Dawn: SPP RICH ~ ACOUSTIC INDICES
mod_df_indices = indices_df %>%
  left_join(y = site_spp_rich, by = 'site') %>%
  drop_na() %>%
  mutate(ARU = factor(substr(site, 4,5)),
         logwavs = log(wavs))

hist(mod_df_indices$site_richness)

mod_df_indices %>%
  select(-site, -logwavs,-ARU, -wavs, -zcr_max, -zcr_min) %>%
  gather(index, value, -site_richness) %>%
  ggplot(aes(x = value, y = site_richness)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~index, scales = "free")

mod1 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Hs, k = 5) + 
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(rugo, k = 5) + 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod1)
par(mfrow = c(2,2))
gam.check(mod1, type = 'response') # Family is okay based on assumptions but upper tail of QQ is low

# drop ARU
mod2 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Hs, k = 5) + 
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(rugo, k = 5) + 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
lrtest(mod1, mod2) # ARU suggested to stay in model
summary(mod2)

# Add ARU/ drop Ht
mod3 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Hs, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(rugo, k = 5) + 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
lrtest(mod1, mod3) # Ht suggested to stay in model
summary(mod3)

# drop ARU/ drop Ht
mod4 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Hs, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(rugo, k = 5) + 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
lrtest(mod1, mod4) # simpler model
summary(mod4)


# drop sfm
mod5 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Hs, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(rugo, k = 5) + 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
lrtest(mod4, mod5) # sfm suggested to stay in
summary(mod5)


# final model
summary(mod4)
par(mfrow = c(2,2))
gam.check(mod4, type = 'response')
draw(mod4, scales = 'fixed', residuals = TRUE)
pred_plot(data_df = mod_df_indices, model_fit = mod4, index_name = 'site_richness')

# check concurvity
# ADI: AEI and Hs
# AEI: ADI and H
# NDSI: NDSI B, NDSI A
# NDSI A: NDSI, NDSI A
# NDSI B: NDSI, NDSI A
# H: Hs and AEI and ZCR
# Hs: H and ADI and ZCR
# ZCR: H and Hs
c <- concurvity(mod4, full = FALSE)$worst

# RMSE
sqrt(mean((mod_df_indices$site_richness - mod4$fitted.values)^2))

sqrt(mean((mod_df_indices$site_richness - mod4$fitted.values)^2))/ 
  diff(range(mod_df_indices$site_richness))

# summarized slopes
ci = slope_summary(mod4)
model_slopes = cbind(index = rep('Bird Spp. Richness', times = nrow(ci)), ci)

# Visualize slope objects 
# calculate all slope values
slope_df = model_slopes

# Use middle 99% of data for each index domain to eliminate spurious tail behavior in plots
slope_df_filtered = slope_df %>%
  select(-smooth) %>%
  group_by(var) %>%
  filter(data >= quantile(data, 0.005),
         data <= quantile(data, 0.995))

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

(gg <- draw(mod4, 
            scales = "fixed", 
            ncol = 4) &
    theme_bw())

#############################################
# Dawn: SPP RICH ~ Biophony + Indices
bio_index_df = mod_df %>%
  left_join(indices_df, by = 'site') %>%
  mutate(Biophony = min_max_norm(Biophony))

bio_index_df %>%
  select(-site, -logwavs, -Anthropophony, -Geophony, -Quiet, -Interference, 
         -sqrtBiophony, -cuberootBiophony, -logBiophony, `wavs.x`, 
         `wavs.y`, -ARU, -zcr_min, -zcr_max) %>%
  gather(index, value, -site_richness) %>%
  ggplot(aes(x = value, y = site_richness)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~index, scales = "free")

mod1 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Hs, k = 5) + 
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(rugo, k = 5) + 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod1)
par(mfrow = c(2,2))
gam.check(mod1, type = 'response')

# drop rugo
mod2 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Hs, k = 5) + 
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
lrtest(mod1, mod2) # rugo suggested to stay
summary(mod1)

# drop Ht
mod3 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Hs, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(rugo, k = 5) + 
             s(sfm, k = 5) + 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
lrtest(mod1, mod3) # Ht stays
summary(mod1)

# drop Hs
mod4 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(rugo, k = 5) + 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
lrtest(mod1, mod4) # drop Hs
summary(mod4)

# drop rugo
mod5 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
lrtest(mod4, mod5) # drop rugo
summary(mod5)

# drop Ht
mod6 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
lrtest(mod5, mod6) # keep Ht
summary(mod5)

# drop ACI
mod7 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
lrtest(mod5, mod7) # keep ACI
summary(mod5)

# drop ARU
mod8 = gam(site_richness ~ 
             logwavs +
             s(Biophony, k = 5) +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
lrtest(mod5, mod8) # drop ARU
summary(mod8)

# drop Ht
mod9 = gam(site_richness ~ 
             logwavs +
             s(Biophony, k = 5) +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
lrtest(mod8, mod9) # keep Ht
summary(mod8)

# drop ACI
mod10 = gam(site_richness ~ 
             logwavs +
             s(Biophony, k = 5) +
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(H, k = 5) + 
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
lrtest(mod8, mod10) # drop ACI
summary(mod10)

# drop ACI
mod11 = gam(site_richness ~ 
              logwavs +
              s(Biophony, k = 5) +
              s(ADI, k = 5) + 
              s(AEI, k = 5) + 
              s(NDSI, k = 5) + 
              s(NDSI_A, k = 5) + 
              s(NDSI_B, k = 5) + 
              s(BI, k = 5) + 
              s(H, k = 5) + 
              s(M, k = 5) + 
              s(R, k = 5) + 
              s(sfm, k = 5) + 
              s(zcr_mean, k = 5),
            data = bio_index_df,
            family = poisson(),
            method = 'ML')
lrtest(mod10, mod11) # keep Ht
summary(mod10)

# drop sfm
mod12 = gam(site_richness ~ 
              logwavs +
              s(Biophony, k = 5) +
              s(ADI, k = 5) + 
              s(AEI, k = 5) + 
              s(NDSI, k = 5) + 
              s(NDSI_A, k = 5) + 
              s(NDSI_B, k = 5) + 
              s(BI, k = 5) + 
              s(H, k = 5) + 
              s(Ht, k = 5) + 
              s(M, k = 5) + 
              s(R, k = 5) + 
              s(zcr_mean, k = 5),
            data = bio_index_df,
            family = poisson(),
            method = 'ML')
lrtest(mod10, mod12) # keep sfm

# drop Ht and sfm
mod13 = gam(site_richness ~ 
              logwavs +
              s(Biophony, k = 5) +
              s(ADI, k = 5) + 
              s(AEI, k = 5) + 
              s(NDSI, k = 5) + 
              s(NDSI_A, k = 5) + 
              s(NDSI_B, k = 5) + 
              s(BI, k = 5) + 
              s(H, k = 5) + 
              s(M, k = 5) + 
              s(R, k = 5) + 
              s(zcr_mean, k = 5),
            data = bio_index_df,
            family = poisson(),
            method = 'ML')
lrtest(mod10, mod13) # drop ht and sfm
summary(mod13)
gam.check(mod13, type = 'response')
draw(mod13, scales = 'fixed', residuals = TRUE)
pred_plot(data_df = bio_index_df, model_fit = mod13, index_name = 'site_richness')

# check concurvity
# ADI: AEI
# AEI: ADI and H
# NDSI: NDSI B, NDSI A
# NDSI A: NDSI, NDSI A
# NDSI B: NDSI, NDSI A
# H: AEI and ZCR
# ZCR: H
c <- concurvity(mod13, full = FALSE)$worst

# RMSE
sqrt(mean((bio_index_df$site_richness - mod13$fitted.values)^2))

sqrt(mean((bio_index_df$site_richness - mod13$fitted.values)^2))/ 
  diff(range(bio_index_df$site_richness))

# summarized slopes
ci = slope_summary(mod13)
model_slopes = cbind(index = rep('Bird Spp. Richness', times = nrow(ci)), ci)

# Visualize slope objects 
# Use middle 95% of data for each index domain to eliminate spurious tail behavior in plots
slope_df_filtered = model_slopes %>%
  select(-smooth) %>%
  group_by(var) %>%
  filter(data >= quantile(data, 0.005),
         data <= quantile(data, 0.995))

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

(gg <- draw(mod13, 
            scales = "fixed", 
            ncol = 3) &
    theme_bw())


ggsave(filename = 'supplementary_figure-Eq5_morning_GAM_PDP.png', 
       plot = gg, 
       device = 'png',
       path = 'figures/',
       width = 8, height = 8, dpi = 500)

# Predicted vs Observed plot
dev <- "Deviance = 65.8%"
bio_index_df$y_pred = predict(mod13, newdata = bio_index_df, type = "response")
(gg <- ggplot(bio_index_df, aes(x = site_richness, y = y_pred)) +
    geom_hex(bins = 25, aes(fill = ..count..)) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    ylab("Predicted Richness") +
    xlab("Observed Richness") +
    theme_bw() +
    scale_fill_viridis_c() +
    index_labeller() +
    theme(text = element_text(size = 16, family = 'Calibri')) +
    annotate("text", x = 10, y = 45, label = dev)
)

ggsave(filename = 'figure_5.png', 
       plot = gg, 
       device = 'png',
       path = 'figures/',
       width = 5, height = 4, dpi = 500)

# Look at potential model saturation
# Compare predicted ~ observed vs log(predicted) ~ observed 
summary(lm(y_pred ~ site_richness, bio_index_df))
summary(lm(log(y_pred) ~ site_richness, bio_index_df))
