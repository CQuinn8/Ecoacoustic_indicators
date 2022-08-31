library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(gratia)
source('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/code/paper1-AcousticIndices/utility_fxs.R')

# Todo:
# 1. rename indices to paper format
# 2. Look at partial effects to see if, for example, rate of change is high for low values of a covariate 
#     so the index is sensitive to ANY presence vs a more steady/linear relationship. E.g. if the curve 
#     is logarithmic any presence of the sound is impactful vs if it is exponential there needs to be a 
#     good amount of the signal to affect the index.

# Formated index names
index_names <- list(
  'ACI' = "ACI",
  'ADI' = "ADI",
  'AEI' = "AEI",
  'BI' = "BI",
  'H' = "H",
  'Hs' = expression(H[s]),
  'Ht' = expression(H[t]),
  'M' = "M",
  'NDSI' = "NDSI",
  'NDSI_A' = "NDSI-\u03B1",
  'NDSI_B' = "NDSI-\u03B2",
  'R' = "R",
  'rugo' = "Rugosity",
  'sfm' = "SFM",
  'zcr_mean' = "ZCR"
)
index_labeller <- function(variable,value){
  return(index_names[value])
}
wd <- '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/'
indices_df <- fread(paste0(wd, 'acoustic_indices_aggregation/averages/site_avg_acoustic_indices.csv'))

# read in abgqi data
abgqi_df <- fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI.csv')) %>%
  mutate(ARU = substr(site, 4,5),
         ARUdevice = substr(site, 4,8)) %>%
  select(-contains('var')) %>%
  select(site, wavs, contains('mean'), ARU, -Unidentified_mean) %>%
  rename(Anthropophony = Anthropophony_mean, Biophony = Biophony_mean, Geophony = Geophony_mean,
         Quiet = Quiet_mean, Interference = Interference_mean)

# Clean sites that are known problems
abgqi_df_no_error <- abgqi_df %>%
  filter(wavs >= (144 * 2.5)) # 44 sites

# sites with static
error_sites <- read.csv('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/data/identified_problem_sites_witherrors.csv')
error_sites <- paste0(error_sites$SiteID, collapse = '|')
mod_df <- abgqi_df_no_error[!grepl(error_sites, abgqi_df_no_error$site),] # 7 sites dropped after filter above

# scale ABGQI
mod_df <- mod_df %>%
  mutate(Anthropophony = min_max_norm(Anthropophony),
         Biophony      = min_max_norm(Biophony),
         Geophony      = min_max_norm(Geophony),
         Quiet         = min_max_norm(Quiet),
         Interference  = min_max_norm(Interference))

# Read in model objects
model_objects <- readRDS(file = paste0(wd, 'modeling/acoustic_indices/LM_model_objects/gam_model_objects_20June2022.RData'))
model_slopes <- readRDS(file = paste0(wd, 'modeling/acoustic_indices/LM_model_objects/gam_model_slopes_20June2022.RData'))



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

length(unique(abgqi_df$site)) # total sites pre-analysis (n = 1,247)
sum(abgqi_df$wavs)            # total wavs pre-analysis (n = 741,061)

mod_df %>%
  select(Anthropophony, Biophony, Geophony, Quiet, Interference) %>%
  dplyr::summarise(across(where(is.numeric), mean_sd)) %>%
  round(6) * 100

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
slope_df_filtered %>%
  group_by(var) %>%
  summarise(min = min(data, na.rm = TRUE),
            max = max(data, na.rm = TRUE))

# Plot slope boxplots
(gg <- slope_df_filtered %>%
    mutate(var = factor(var, levels = c( 'Interference','Quiet', 'Geophony','Biophony','Anthropophony'))) %>%
    ggplot(aes(x = factor(var), y = derivative)) +
    geom_boxplot(width = 0.5, 
                 outlier.shape = NA, 
                 position = position_dodge(preserve = "single"),
                 notch = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ylab("Slope") +
    xlab("Soundscape Source") +
    scale_fill_manual(values = c('#B8B8B8', '#FFFFFF'), name = "") +
    coord_flip() +
    theme_bw() +
    facet_wrap(~index, scales = 'free_x', labeller = index_labeller) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    theme(legend.position = 'bottom',
          text = element_text(size = 16, family = 'Calibri'),
          axis.text.x = element_text(angle = 45, hjust = 1)))

# (gg = annotate_figure(gg, top = text_grob("95% CI for first derivative (slopes) of GAM partial effects.
# Note: should be interpreted alongside partial effects plots.", size = 14, face = "bold"),
#                 bottom = text_grob("Values reflect middle 99% of covariate range to minimize erroneous tail behavior.",
#                                    size = 10)))

ggsave(gg, filename = 'G:/My Drive/NAU/Dissertation/paper2-AcousticIndices/figures/FigR1-GAM_slopes_by_model.png', 
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

# Table R1
R1 = Reduce(function(...) merge(..., all = TRUE, by = 'variable'), 
       list(formulas, links, rsq, deviance, aics))
R1 = R1[order(R1$Deviance, decreasing = TRUE),]
print(R1$variable)
# Organize table for output
library(kableExtra)
R1 %>%
  kbl(caption = "Table R1") %>%
  row_spec(0, bold = T) %>%
  kable_classic(full_width = F, html_font = "Times New Roman")
  
  #knitr::kable(digits = 4, caption = "Results Table R1")


######################################
# Visualize partial effects + interpretation
temp_df = merge(x = indices_df, y = mod_df)
  
# y-axis is the centered smooth
draw(model_objects$zcr_mean, scales = "fixed", residuals = TRUE)
# most affected by anthro, negatively
# Bio has small effect mostly at lower values too
# Geophony and int are negatively related
# Quiet is mixed

draw(model_objects$ACI, scales = "fixed", residuals = TRUE)
model_objects$ACI
# most affected by interference esp lower values - makes sense as these are staccato 
# Biophony has a positively linear relationship
# Geophony also is positive (slightly more than Bio)
# Anthro is increasing at lesser rate

draw(model_objects$M, scales = "fixed", residuals = TRUE)
# minimal biophony effect
# strong log shaped Anthro and exponential Interference
# Geophony is positive until some mixed behavior at sparser mid-high values

draw(model_objects$H, scales = "fixed", residuals = TRUE)
summary(model_objects$H)

draw(model_objects$Hs, scales = "fixed", residuals = TRUE)
summary(model_objects$Hs)

draw(model_objects$NDSI_A, scales = "fixed", residuals = TRUE)
summary(model_objects$NDSI_A)

draw(model_objects$AEI, scales = "fixed", residuals = TRUE)

draw(model_objects$Ht, scales = "fixed", residuals = TRUE)
summary(model_objects$Ht)

draw(model_objects$rugo, scales = "fixed", residuals = TRUE)

draw(model_objects$NDSI_B, scales = "fixed", residuals = TRUE)
summary(model_objects$NDSI_B)

draw(model_objects$NDSI, scales = "fixed", residuals = TRUE)
summary(model_objects$NDSI)
# weakly impacted by interference (positive)
# negative Anthro
# Strong positive effect from Bio and lessser so Quiet

draw(model_objects$ADI, scales = "fixed", residuals = TRUE)

draw(model_objects$sfm, scales = "fixed", residuals = TRUE)
model_objects$sfm
summary(indices_df$sfm)

draw(model_objects$BI, scales = "fixed", residuals = TRUE)
pred_plot(data_df = temp_df, model_fit = model_objects$BI, index_name = "BI")

draw(model_objects$R, scales = "fixed", residuals = TRUE)
