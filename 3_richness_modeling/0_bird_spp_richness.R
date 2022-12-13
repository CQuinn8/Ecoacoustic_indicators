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
wd = 'zenodo/'

# bird species pres/abs by site 
bird_df = fread('data/sonoma_s2l_predictions_ROImaxF05_220905_summarized.csv')

# reduce sites used in ABGQI analyses
sites = fread("data/sites_used_in_GAMs.csv")

# ABGQI
abgqi_df = fread('data/site_avg_ABGQI.csv') %>%
  mutate(ARU = substr(site, 4,5),
         ARUdevice = substr(site, 4,8)) %>%
  select(-contains('var')) %>%
  select(site, wavs, contains('mean'), ARU, -Unidentified_mean) %>%
  rename(Anthropophony = Anthropophony_mean, Biophony = Biophony_mean, Geophony = Geophony_mean,
         Quiet = Quiet_mean, Interference = Interference_mean)

# Acoustic Indices
indices_df = fread('data/site_avg_acoustic_indices.csv')

# sites with static
error_sites = read.csv('data/identified_problem_sites_witherrors.csv')
error_sites = paste0(error_sites$SiteID, collapse = '|')

###################################
# Prep data and summarize
# subset ABGQI sites results in 1226 -> 1185 sites
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

# Clean sites that are known problems
abgqi_df = abgqi_df %>%
  filter(wavs >= (144 * 2.5)) # 44 sites
abgqi_df_no_error = abgqi_df[!grepl(error_sites, abgqi_df$site),]

# subset indices to target sites
indices_df = indices_df[indices_df$site %in% sites$sites, ]

# Merge dataframes dropping the 10 sites with no bird spp observations
mod_df = abgqi_df_no_error %>%
  left_join(y = site_spp_rich, by = 'site') %>%
  drop_na()
  #replace_na(list(site_richness = 0))

mean(mod_df$site_richness)
sd(mod_df$site_richness)

mod_df = mod_df %>%
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) 

mean(mod_df$Biophony)# 0.2398
sd(mod_df$Biophony)  # 0.1546

mod_df %>%
  group_by(ARU) %>%
  summarise(n())

# rate of bird presences 
bird_df_nwavs <- abgqi_df_no_error %>%
  select(site, wavs) %>%
  left_join(bird_df, by = 'site') %>%
  drop_na()

rate_birds <- bird_df_nwavs %>%
  mutate_at(vars(-site,-wavs), funs(. / wavs)) %>%
  summarise(across(where(is.numeric), mean))

# write.csv(rate_birds, row.names = F, paste0(wd,'results/modeling/bird_rate_24hr.csv'))

#############################################
# 24hr: BIOPHONY ~ SPP RICH 
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
# 24hr: SPP RICH ~ ACOUSTIC INDICES
mod_df_indices = indices_df %>%
  left_join(y = site_spp_rich, by = 'site') %>%
  drop_na()
mod_df_indices = merge(mod_df_indices, abgqi_df_no_error, by = 'site') %>%
  mutate(logwavs = log(wavs),
         ARU = factor(ARU),
         logRichness = log(site_richness+1e-7)) %>%
  select(-YYYY, -site, -wavs, -Anthropophony, -Biophony, -Geophony, -Quiet, -Interference)

hist(mod_df_indices$site_richness)
hist(mod_df_indices$logRichness)

mod_df_indices %>%
  select(-logwavs, -logRichness) %>%
  gather(index, value, -ARU, -site_richness) %>%
  ggplot(aes(x = value, y = site_richness, colour = ARU)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam") +
  facet_wrap(~index, scales = "free")

# extreme low outlier spp richness = 1
mod_df_indices = mod_df_indices %>%
  filter(site_richness > 0)

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

# drop ARU
mod2 = gam(site_richness ~ 
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
             s(rugo, k = 5)+ 
             s(sfm, k = 5)+ 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod2)
lrtest(mod1, mod2) # ARU not recommended to drop

# drop rugo/add ARU
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
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod3)
lrtest(mod2, mod3)

# Test ARU again
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
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) + 
             s(sfm, k = 5) + 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
lrtest(mod3, mod4) # keep ARU

# drop ADI
mod5 = gam(site_richness ~ 
             logwavs +
             ARU +
             s(ACI, k = 5) + 
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
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod5)
lrtest(mod3, mod5)

# test ARU
mod6 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
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
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
lrtest(mod5, mod6) # drop
summary(mod6)

# drop H 
mod7 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) + 
             s(NDSI_A, k = 5) + 
             s(NDSI_B, k = 5) + 
             s(BI, k = 5) + 
             s(Hs, k = 5) + 
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) +
             s(sfm, k = 5) + 
             s(zcr_mean, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod7)
lrtest(mod6, mod7)

# drop zcr
mod8 = gam(site_richness ~ 
             logwavs +
             s(ACI, k = 5) + 
             s(AEI, k = 5) + 
             s(NDSI, k = 5) +
             s(NDSI_A, k = 5) +
             s(NDSI_B, k = 5) +
             s(BI, k = 5) + 
             s(Hs, k = 5) + 
             s(Ht, k = 5) + 
             s(M, k = 5) + 
             s(R, k = 5) +
             s(sfm, k = 5),
           data = mod_df_indices,
           family = poisson(),
           method = 'ML')
summary(mod8)
lrtest(mod7, mod8) # keep zcr - mod7 is final

# final model
summary(mod7)
par(mfrow = c(2,2))
gam.check(mod7, type = 'response')
draw(mod7, scales = 'fixed', residuals = TRUE)
pred_plot(data_df = mod_df_indices, model_fit = mod7, index_name = 'site_richness')

# check concurvity - NDSI, NDSI_A, NDSI_B
# NDSI only: sign changes from negative to positive
# NDSI_B only: remains relatively stable
# NDSI_A only: remains relatively stable
concurvity(mod7, full = TRUE)
c <- concurvity(mod7, full = FALSE)$worst 

# RMSE
sqrt(mean((mod_df_indices$site_richness - mod7$fitted.values)^2))
sqrt(mean((mod_df_indices$site_richness - mod7$fitted.values)^2))/ 
  diff(range(mod_df_indices$site_richness))

# parametric wavs effect
exp(0.2988)

# summarized slopes
ci = slope_summary(mod7)
model_slopes = cbind(index = rep('Bird Spp. Richness', times = nrow(ci)), ci)

# Visualize slope objects 
# calculate all slope values
slope_df = model_slopes

# Use middle 95% of data for each index domain
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

(gg <- draw(mod6, 
            scales = "fixed", 
            ncol = 4) &
    theme_bw())


#############################################
# 24hr: SPP RICH ~ Biophony
mod_df_bio = mod_df %>%
  mutate(ARU = factor(ARU),
         logwavs = log(wavs)) %>%
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

sqrt(mean((mod_df_bio$site_richness - mod1$fitted.values)^2))/ 
  diff(range(mod_df_bio$site_richness))

mod_df$spp_rich_pred = predict(mod1, newdata = mod_df_bio, type = "response")
pred_plot(mod_df, mod1, "site_richness")

#############################################
# 24hr: SPP RICH ~ Biophony + Indices
bio_index_df = mod_df %>%
  left_join(indices_df, by = 'site') %>%
  mutate(Biophony = min_max_norm(Biophony))

bio_index_df %>%
  select(-site, -logwavs, -Anthropophony, -Geophony, -Quiet, -Interference, 
         -sqrtBiophony, -cuberootBiophony, -logBiophony, -YYYY) %>%
  gather(index, value, -ARU, -site_richness) %>%
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

# drop rugo
mod2 = gam(site_richness ~ 
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
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod2)
gam.check(mod2, type = 'response')
lrtest(mod1, mod2)

# drop sfm 
mod3 = gam(site_richness ~ 
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
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod3)
gam.check(mod3, type = 'response')
lrtest(mod2, mod3)

# drop ht 
mod4 = gam(site_richness ~ 
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
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod4)
gam.check(mod4, type = 'response')
lrtest(mod3, mod4) # not suggested to drop Ht
summary(mod3)

# add Ht, drop hs
mod5 = gam(site_richness ~ 
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
             s(Ht, k = 5)+
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod5)
gam.check(mod5, type = 'response')
lrtest(mod3, mod5) # Hs recommended to drop

# drop Ht 
mod6 = gam(site_richness ~ 
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
             s(M, k = 5)+ 
             s(R, k = 5)+ 
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod6)
lrtest(mod5, mod6) # keep Ht
summary(mod5)

# test ARU 
mod7 = gam(site_richness ~ 
             logwavs +
             s(Biophony, k = 5) +
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
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod7)
lrtest(mod5, mod7) # keep ARU
summary(mod5)

# test R 
mod8 = gam(site_richness ~ 
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
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod8)
lrtest(mod5, mod8)

# drop Ht
mod9 = gam(site_richness ~ 
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
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod9)
lrtest(mod8, mod9)

# drop ARU
mod10 = gam(site_richness ~ 
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
             s(zcr_mean, k = 5),
           data = bio_index_df,
           family = poisson(),
           method = 'ML')
summary(mod10)
lrtest(mod9, mod10) # suggest to retain ARU
summary(mod9) # final model

# Final model - mod9
summary(mod9)
draw(mod9, scales = 'fixed', residuals = TRUE)
pred_plot(data_df = bio_index_df, model_fit = mod5, index_name = 'site_richness')

# check concurvity
# AEI, H, and ADI; NDSI, NDSI_A, NDSI_B: ZCR and H
concurvity(mod9, full = TRUE)
c <- concurvity(mod9, full = FALSE)$worst

# RMSE
sqrt(mean((bio_index_df$site_richness - mod9$fitted.values)^2))

sqrt(mean((bio_index_df$site_richness - mod9$fitted.values)^2))/ 
  diff(range(bio_index_df$site_richness))

# summarized slopes
ci = slope_summary(mod9)
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

(gg <- draw(mod9, 
            scales = "fixed", 
            ncol = 3) &
    theme_bw())

# Investigate why at low richness we overpredict compared to higher (possibly due to gradient in disturbance)
mod_df %>%
  select(Anthropophony, Quiet, site_richness) %>%
  pivot_longer(names_to = "sound", values_to = "value", !site_richness) %>%
  ggplot(aes(x = value, y = site_richness)) +
    geom_hex() +
    facet_wrap(. ~ sound)
