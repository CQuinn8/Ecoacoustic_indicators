library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(gratia)
library(kableExtra)
source('utility_fxs.R')

wd <- '/paper-AcousticIndices/results/'
out_dir <- '/figures/'

# Site average acoustic indices
indices_df <- fread(paste0(wd, 'acoustic_indices_aggregation/averages/site_avg_acoustic_indices.csv'))

# Site Avg ABGQI data
abgqi_df <- fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI.csv')) %>%
  mutate(ARU = substr(site, 4,5),
         ARUdevice = substr(site, 4,8)) %>%
  select(-contains('var')) %>%
  select(site, wavs, contains('mean'), ARU, -Unidentified_mean) %>%
  rename(Anthropophony = Anthropophony_mean, Biophony = Biophony_mean, Geophony = Geophony_mean,
         Quiet = Quiet_mean, Interference = Interference_mean)

# Clean sites that are known problems
abgqi_df_no_error <- abgqi_df %>%
  filter(wavs >= (144 * 2.5)) # rm 44 sites

# sites with static
error_sites <- read.csv('/paper-AcousticIndices/data/identified_problem_sites_witherrors.csv')
error_sites <- paste0(error_sites$SiteID, collapse = '|')
mod_df <- abgqi_df_no_error[!grepl(error_sites, abgqi_df_no_error$site),] # 7 sites dropped after filter above

# scale ABGQI to match models
mod_df <- mod_df %>%
  mutate(Anthropophony = min_max_norm(Anthropophony),
         Biophony      = min_max_norm(Biophony),
         Geophony      = min_max_norm(Geophony),
         Quiet         = min_max_norm(Quiet),
         Interference  = min_max_norm(Interference))

# Read in model objects
model_objects <- readRDS(file = '/models/gam_model_objects.RData')
model_slopes <- readRDS(file = '/models/gam_model_slopes.RData')

######################################
# Data summaries
length(unique(mod_df$site)) # sites used in analyses
sum(mod_df$wavs)            # number of recordings in analyses

# count number of sites per year (note - one site mislabeled 2041 instead of 2021)
mod_df %>%
  mutate(YY = substr(site, 10, 11)) %>%
  group_by(YY) %>%
  count()

# get number of recordings per year
mod_df %>%
  mutate(YY = substr(site, 10, 11)) %>%
  group_by(YY) %>%
  summarise(sum(wavs))

# Pre-analysis sites and recordings
length(unique(abgqi_df$site)) # total sites pre-analysis (n = 1,247)
sum(abgqi_df$wavs)            # total wavs pre-analysis (n = 741,061)

# Avg percentages for ABGQI
mod_df %>%
  select(Anthropophony, Biophony, Geophony, Quiet, Interference) %>%
  dplyr::summarise(across(where(is.numeric), mean_sd)) %>%
  round(6) * 100

# Visual number of sites by year
mod_df %>%
  mutate(YY = substr(site, 10, 11)) %>%
  filter(YY != 41) %>%
  ggplot(aes(x = YY)) +
    geom_bar(fill = '#3E97B6') +
    theme_bw() +
    ggtitle('Number of survey sites per year') +
    xlab('Year') +
    ylab('Count of sites')
  

# Visual number of recordings by year
mod_df %>%
  mutate(YY = substr(site, 10, 11)) %>%
  filter(YY != 41) %>%
  group_by(YY) %>%
  mutate(Count = sum(wavs) / 726378 * 1195) %>%
  ggplot(aes(x = YY)) +
    geom_bar(fill = '#3E97B6') +
  geom_point(aes(y = Count, x = YY), size = 5) +
    theme_bw() +
    ggtitle('Number of survey sites per year & 
            proportion of annual recordings') +
    xlab('Year') +
    ylab('Count of sites')

######################################
# Visualize slope objects
# calculate all slope values
slope_df <- do.call("rbind", model_slopes)

# use middle 99% of x values to plot slopes to avoid extreme tail behavior
abgqi_limits <- mod_df %>%
  select(Anthropophony, Biophony, Geophony, Quiet, Interference) %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(lower = quantile(value, 0.005),
            upper = quantile(value, 0.995))

# filter slope df based on percentile values
slope_df_filtered <- slope_df %>%
  filter(case_when(var == 'Anthropophony' ~ data >= abgqi_limits$lower[1] & data <= abgqi_limits$upper[1],
                   var == 'Biophony'      ~ data >= abgqi_limits$lower[2] & data <= abgqi_limits$upper[2],
                   var == 'Geophony'      ~ data >= abgqi_limits$lower[3] & data <= abgqi_limits$upper[3],
                   var == 'Interference'  ~ data >= abgqi_limits$lower[4] & data <= abgqi_limits$upper[4],
                   var == 'Quiet'         ~ data >= abgqi_limits$lower[5] & data <= abgqi_limits$upper[5]))

# Summarize the filtered ABGQI values
slope_df_filtered %>%
  group_by(var) %>%
  summarise(min = min(data, na.rm = TRUE),
            max = max(data, na.rm = TRUE))

# Plot slope boxplots (i.e., Fig 2)
(gg <- slope_df_filtered %>%
    # relevel factors
    mutate(var = factor(var, levels = c( 'Interference','Quiet', 'Geophony','Biophony','Anthropophony'))) %>%
    
    # Instantiate ggplot
    ggplot(aes(x = factor(var), y = derivative)) +
      geom_boxplot(width = 0.5, 
                   outlier.shape = NA, 
                   position = position_dodge(preserve = "single"),
                   notch = TRUE) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ylab("Slope") +
      xlab("Soundscape Components") +
      scale_fill_manual(values = c('#B8B8B8', '#FFFFFF'), name = "") +
      coord_flip() +
      theme_bw() +
    
      # use label key in utility_fxs.R to correctly label Acoustic Indices
      facet_wrap(~index, scales = 'free_x', labeller = index_labeller) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
      theme(legend.position = 'bottom',
            text = element_text(size = 16, family = 'Calibri'),
            axis.text.x = element_text(angle = 45, hjust = 1)))

# Save figure 2
ggsave(gg, filename = paste0(out_dir,'figure2.png'), 
       height = 6, width = 6.5, unit = 'in', dpi = 500)

######################################
# Summarize model fits
temp = summary(model_objects[[1]])

# adj R-sq values
rsq = data.frame(do.call(cbind, lapply(model_objects, function(x) summary(x)$r.sq))) %>%
  gather(variable, Rsq)
mean(rsq$Rsq)
sd(rsq$Rsq)

# deviance values
deviance = data.frame(do.call(cbind, lapply(model_objects, function(x) summary(x)$dev.expl))) %>%
  gather(variable, Deviance)
mean(deviance$Deviance)
sd(deviance$Deviance)

# model family
formulas = data.frame(do.call(cbind, lapply(model_objects, function(x) summary(x)$family$family))) %>%
  gather(variable, family)

# AIC 
aics = data.frame(do.call(cbind, lapply(model_objects, function(x) AIC(x)))) %>%
  gather(variable, AIC)

# links
links = data.frame(do.call(cbind, lapply(model_objects, function(x) summary(x)$family$link))) %>%
  gather(variable, link)

# Table 3
# Combine the separately created numeric summaries
R1 = Reduce(function(...) merge(..., all = TRUE, by = 'variable'), 
       list(formulas, links, rsq, deviance, aics))

# Order the combined metrics by decreasing deviance
R1 = R1[order(R1$Deviance, decreasing = TRUE),]
print(R1$variable)

# Organize table for output
R1 %>%
  kbl(caption = "Table 3") %>%
  row_spec(0, bold = T) %>%
  kable_classic(full_width = F, html_font = "Times New Roman")


######################################
# Parital Dependence Plots (Supplementary Material)
# Each acoustic index PDP summary is drawn and then saved using the write_pdp() in utility_fxs.R
#  - Some indices have summary notes included
#  - Oganized in descending deviance

# combine dataframes
temp_df = merge(x = indices_df, y = mod_df)

# ZCR
(gg <- draw(model_objects$zcr_mean, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
write_pdp(plot_obj = gg, index_name = "zcr", rows = 2, out_dir = out_dir)
summary(model_objects$zcr_mean)
# most affected by anthro, negatively
# Bio has small effect mostly at lower values too
# Geophony and int are negatively related
# Quiet is mixed

# ACI
(gg <- draw(model_objects$ACI, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
write_pdp(plot_obj = gg, index_name = "aci", rows = 1, out_dir = out_dir)
model_objects$ACI
# most affected by interference esp lower values - makes sense as these are staccato 
# Biophony has a positively linear relationship
# Geophony also is positive (slightly more than Bio)
# Anthro is increasing at lesser rate

# M
(gg <- draw(model_objects$M, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
write_pdp(plot_obj = gg, index_name = "M", rows = 2, out_dir = out_dir)
# minimal biophony effect
# strong log shaped Anthro and exponential Interference
# Geophony is positive until some mixed behavior at sparser mid-high values

# H
(gg <- draw(model_objects$H, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
write_pdp(plot_obj = gg, index_name = "H", rows = 2, out_dir = out_dir)
summary(model_objects$H)

# Hs
(gg <- draw(model_objects$Hs, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
write_pdp(plot_obj = gg, index_name = "Hs", rows = 2, out_dir = out_dir)
summary(model_objects$Hs)

# NDSI_A
(gg <- draw(model_objects$NDSI_A, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
write_pdp(plot_obj = gg, index_name = "NDSI_A", rows = 1, out_dir = out_dir)
summary(model_objects$NDSI_A)

# AEI
(gg <- draw(model_objects$AEI, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
summary(model_objects$AEI)
write_pdp(plot_obj = gg, index_name = "AEI", rows = 2, out_dir = out_dir)

# Ht
(gg <- draw(model_objects$Ht, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
write_pdp(plot_obj = gg, index_name = "Ht", rows = 2, out_dir = out_dir)
summary(model_objects$Ht)

# Rugosity
(gg <- draw(model_objects$rugo, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
summary(model_objects$rugo)
write_pdp(plot_obj = gg, index_name = "rugo", rows = 2, out_dir = out_dir)

# NDSI_B
(gg <- draw(model_objects$NDSI_B, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
summary(model_objects$NDSI_B)
write_pdp(plot_obj = gg, index_name = "NDSI_B", rows = 1, out_dir = out_dir)

# NDSI
(gg <- draw(model_objects$NDSI, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
summary(model_objects$NDSI)
write_pdp(plot_obj = gg, index_name = "NDSI", rows = 1, out_dir = out_dir)
# weakly impacted by interference (positive)
# negative Anthro
# Strong positive effect from Bio and lessser so Quiet

# ADI
(gg <- draw(model_objects$ADI, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
summary(model_objects$ADI)
write_pdp(plot_obj = gg, index_name = "ADI", rows = 2, out_dir = out_dir)

# SFM
(gg <- draw(model_objects$sfm, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
summary(model_objects$sfm)
write_pdp(plot_obj = gg, index_name = "SFM", rows = 2, out_dir = out_dir)
model_objects$sfm

# BI
(gg <- draw(model_objects$BI, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
summary(model_objects$BI)
write_pdp(plot_obj = gg, index_name = "BI", rows = 1, out_dir = out_dir)
# pred_plot(data_df = temp_df, model_fit = model_objects$BI, index_name = "BI")

# R
(gg <- draw(model_objects$R, scales = "fixed", residuals = TRUE, ncol = 4) & theme_bw())
summary(model_objects$R)
write_pdp(plot_obj = gg, index_name = "R", rows = 1, out_dir = out_dir)