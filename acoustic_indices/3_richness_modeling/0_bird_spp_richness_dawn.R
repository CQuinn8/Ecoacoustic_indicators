
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

###################################
# DATA IMPORT
wd = '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/'

# bird species pres/abs by site
bird_df = fread('G:/Shared drives/NASA  S2L/Paper Development/Sonoma SDM update/bird_cnn_predictions/sonoma_s2l_predictions_ROImaxF05_220905_summarized-dawn_4am-12pm.csv')

# reduce sites used in ABGQI analyses
sites = fread("//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/sites_used_in_GAMs_21July22.csv")

# ABGQI
abgqi_df = fread(paste0(wd, 'results/ABGQI_inference/averages/site_ABGQI_dawn_4am-12pm.csv')) %>%
  mutate(ARU = substr(site, 4,5)) %>%
  select(-Unidentified)

# Acoustic Indices
indices_df = fread(paste0(wd, 'results/acoustic_indices_aggregation/averages/site_acoustic_indices_dawn_4am-12pm.csv'))

# sites with static
error_sites = read.csv('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/data/identified_problem_sites_witherrors.csv')
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

mean(mod_df$site_richness)
sd(mod_df$site_richness)

mod_df = mod_df %>%
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) 

# rate of bird presences 
bird_df_nwavs <- abgqi_df_no_error %>%
  select(site, wavs) %>%
  left_join(bird_df, by = 'site') %>%
  drop_na()

rate_birds <- bird_df_nwavs %>%
  mutate_at(vars(-site,-wavs), funs(. / wavs)) %>%
  summarise(across(where(is.numeric), mean))

write.csv(rate_birds, row.names = F, paste0(wd,'results/modeling/bird_rate_dawn.csv'))

bird_rate_24hr <- read.csv(paste0(wd,'results/modeling/bird_rate_24hr.csv'))
bird_rate_24hr$dataset <- "24hr"
rate_birds$dataset <- "dawn"

combined_bird_rate <- bird_rate_24hr %>%
  bind_rows(rate_birds)

write.csv(combined_bird_rate, row.names = F, paste0(wd,'results/modeling/bird_rate_combined.csv'))

#############################################
# Dawn: BIOPHONY ~ SPP RICH 
# visualize data
mod_df %>%
  select(Biophony, site_richness, ARU) %>%
  ggplot(aes(y = Biophony, x = site_richness, colour = ARU)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "gam")
mod_df$sqrtBiophony = sqrt(mod_df$Biophony)
mod_df$cuberootBiophony = (mod_df$Biophony)^(1/3)
mod_df$logBiophony = log(mod_df$Biophony+1e-7)

hist(mod_df$Biophony)
hist(mod_df$sqrtBiophony)
hist(mod_df$cuberootBiophony)
hist(mod_df$logBiophony)
# mod_df = mod_df[-886,]

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

# RMSE
sqrt(mean((mod_df$Biophony - modBio1$fitted.values^3)^2))

# Parameteric values
# ARUlg: -0.03512
fromLogit(-0.03512) # LG ARUs have smaller effect on Biophony than AM

# logwavs -0.04519
fromLogit(-0.04519) # logwavs, as number of wavs increases Biophony decreases

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

ggsave(filename = 'figure_r4.png', 
       plot = gg, 
       device = 'png',
       path = 'G:/My Drive/NAU/Dissertation/paper2-AcousticIndices/figures/',
       width = 6, height = 4, dpi = 500)

# # grouped by LULC
# mod_df_remod %>%
#   select(site_richness, Biophony, Biophony_pred, ARU, LULC) %>%
#   gather(covariate, value, -Biophony, -Biophony_pred, -ARU, -LULC) %>%
#   ggplot(aes(x = value, y = Biophony_pred)) +
#   geom_point(aes(x = value, y = Biophony, colour = LULC), alpha = 0.3, size = 2,
#              position = position_jitter(w = 0.1, h = 0)) +
#   #geom_point(alpha = 0.1) +
#   # geom_smooth(method = 'gam', colour = 'black', alpha = 0.4) +
#   geom_smooth(method = 'glm', aes(colour = LULC), alpha = 0.4, se = FALSE) +
#   ggtitle("Biophony: black = predicted values") +
#   labs(y = 'Percent Biophony', x = 'Bird Species Richness') +
#   theme_bw() +
#   scale_colour_brewer(palette = "Set1")
# 
# pred_plot(data_df = mod_df_remod, model_fit = remodBio1, index_name = "cuberootBiophony")
# 
# mod_df_remod %>%
#   select(LULC, site_richness) %>%
#   group_by(LULC) %>%
#   summarise(mean(site_richness),
#             sd(site_richness),
#             n(site_richness))
# pred_plot(data_df = mod_df_remod, model_fit = remodBio1, index_name = "cuberootBiophony")

#############################################
# Dawn: SPP RICH ~ ACOUSTIC INDICES
mod_df_indices = indices_df %>%
  left_join(y = site_spp_rich, by = 'site') %>%
  drop_na() %>%
  mutate(ARU = factor(substr(site, 4,5)),
         logwavs = log(wavs))

hist(mod_df_indices$site_richness)

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
par(mfrow = c(2,2))
gam.check(mod1, type = 'response') # Family is okay based on assumptions but upper tail of QQ is too low

# drop ADI (p = 0.529434)
mod2 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(ACI, k = 5) + 
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
summary(mod2)
AIC(mod1, mod2)

# drop H (p = 0.47638)
mod3 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(ACI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
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
summary(mod3)
AIC(mod2, mod3)

# drop ARU (p = 0.305)
mod4 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
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
summary(mod4)
AIC(mod3, mod4)

# drop sfm (p = 0.325077)
mod5 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(Hs, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(rugo, k = 5)+ 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod5)
AIC(mod4, mod5)

# drop Rugo (p = 0.16588)
mod6 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(Hs, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod6)
AIC(mod5, mod6) # dAIC < 6 so more parsimonious model suggested

# test Hs (p = 0.07418)
mod7 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod7)
AIC(mod6, mod7) # dAIC < 6 so more parsimonious model suggested
AIC(mod5, mod7) # dAIC < 6 so more parsimonious model suggested

# test zcr (p = 0.322121)
mod8 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod8)
AIC(mod7, mod8) # dAIC < 6 so more parsimonious model suggested
AIC(mod6, mod8) # dAIC < 6 so more parsimonious model suggested
AIC(mod5, mod8) # dAIC < 6 so more parsimonious model suggested

# test Ht (p = 0.182424)
mod9 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod9)
AIC(mod8, mod9) # dAIC < 6 so more parsimonious model suggested
AIC(mod5, mod9) # dAIC < 6 so more parsimonious model suggested

# mod
summary(mod9)
par(mfrow = c(2,2))
gam.check(mod9, type = 'response')
draw(mod7, scales = 'fixed', residuals = TRUE)
pred_plot(data_df = mod_df_indices, model_fit = mod9, index_name = 'site_richness')

# RMSE
sqrt(mean((mod_df_indices$site_richness - mod9$fitted.values)^2))

sqrt(mean((mod_df_indices$site_richness - mod9$fitted.values)^2))/ 
  diff(range(mod_df_indices$site_richness))

# parametric wavs effect
exp(0.2988)

# summarized slopes
ci = slope_summary(mod6)
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

(gg <- draw(mod6, 
            scales = "fixed", 
            ncol = 4) &
    theme_bw())

ggsave(filename = 'figure_6.png', 
       plot = gg, 
       device = 'png',
       path = 'G:/My Drive/NAU/Dissertation/paper2-AcousticIndices/figures/',
       width = 8, height = 8, dpi = 500)

#############################################
# DAwn: SPP RICH ~ Biophony
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

# RMSE
sqrt(mean((mod_df_bio$site_richness - mod1$fitted.values)^2))

# rRMSE 
sqrt(mean((mod_df_bio$site_richness - mod1$fitted.values)^2)) / 
       diff(range(mod_df_bio$site_richness))


mod_df_bio$spp_rich_pred = predict(mod1, newdata = mod_df_bio, type = "response")
pred_plot(mod_df_bio, mod1, "site_richness")

#############################################
# Dawn: SPP RICH ~ Biophony + Indices
bio_index_df = mod_df %>%
  left_join(indices_df, by = 'site') %>%
  mutate(Biophony = min_max_norm(Biophony))

mod1 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
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
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod1)
par(mfrow = c(2,2))
gam.check(mod1, type = 'response')

# drop NDSI (p = 0.879856)
mod2 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ACI, k = 5) + 
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
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
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod2)
AIC(mod1, mod2)

# drop ACI (p = 0.821472)
mod3 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
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
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod3)
AIC(mod2, mod3)


# drop Hs (p = 0.764678)
mod4 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(sfm, k = 5)+ 
             s(rugo, k = 5)+ 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod4)
AIC(mod3, mod4)

# drop Rugo (p = 0.428582)
mod5 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(sfm, k = 5)+ 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod5)
AIC(mod4, mod5)

# drop SFM (p = 0.471827)
mod6 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(Biophony, k = 5) +
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod6)
AIC(mod5, mod6)

# drop ARU (p = 0.215 )
mod7 = gam(site_richness ~ 
             logwavs +
             s(Biophony, k = 5) +
             s(ADI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod7)
AIC(mod6, mod7)

# drop ADI (p = 0.404252)
mod8 = gam(site_richness ~ 
             logwavs +
             s(Biophony, k = 5) +
             s(AEI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod8)
AIC(mod7, mod8)

# drop R (p = 0.197002)
mod9 = gam(site_richness ~ 
             logwavs +
             s(Biophony, k = 5) +
             s(AEI, k = 5) + 
             s(NDSI_A, k = 5)+ 
             s(NDSI_B, k = 5)+ 
             s(BI, k = 5)+ 
             s(H, k = 5)+ 
             s(Ht, k = 5)+ 
             s(M, k = 5)+ 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod9)
AIC(mod8, mod9)


# test Ht (p = 0.056907)
mod10 = gam(site_richness ~ 
              logwavs +
              s(Biophony, k = 5) +
              s(AEI, k = 5) + 
              s(NDSI_A, k = 5)+ 
              s(NDSI_B, k = 5)+ 
              s(BI, k = 5)+ 
              s(H, k = 5)+ 
              s(M, k = 5)+ 
              s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod10)
AIC(mod9, mod10) # dAIC < 6 - use more parsimonious model 10
gam.check(mod10, type = 'response')

draw(mod10, scales = 'fixed', residuals = TRUE)
pred_plot(data_df = bio_index_df, model_fit = mod10, index_name = 'site_richness')

# RMSE
sqrt(mean((bio_index_df$site_richness - mod10$fitted.values)^2))

sqrt(mean((bio_index_df$site_richness - mod10$fitted.values)^2))/ 
  diff(range(bio_index_df$site_richness))

# summarized slopes
ci = slope_summary(mod10)
model_slopes = cbind(index = rep('Bird Spp. Richness', times = nrow(ci)), ci)

# Visualize slope objects 
# Use middle 95% of data for each index domain to eliminate spurious tail behavior in plots
slope_df_filtered = model_slopes %>%
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

(gg <- draw(mod10, 
            scales = "fixed", 
            ncol = 3) &
    theme_bw())

ggsave(filename = 'figure_6.png', 
       plot = gg, 
       device = 'png',
       path = 'G:/My Drive/NAU/Dissertation/paper2-AcousticIndices/figures/',
       width = 8, height = 8, dpi = 500)

# Predicted vs Observed plot
bio_index_df$y_pred = predict(mod10, newdata = bio_index_df, type = "response")
(gg = ggplot(bio_index_df, aes(x = site_richness, y = y_pred)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    ylab("Predicted Richness") +
    xlab("Observed Richness") +
    theme_bw() +
    theme(text = element_text(size = 16, family = 'Calibri')))


ggsave(filename = 'figure_5.png', 
       plot = gg, 
       device = 'png',
       path = 'G:/My Drive/NAU/Dissertation/paper2-AcousticIndices/figures/',
       width = 4, height = 4, dpi = 500)
