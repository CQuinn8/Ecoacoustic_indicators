
# Interpret deviance and Rsq: https://stats.stackexchange.com/questions/190172/how-i-can-interpret-gam-results


library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
source('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/code/paper1-AcousticIndices/utility_fxs.R')

# Todo:
# 1. rename indices to paper format
# 2. Look at partial effects to see if, for example, rate of change is high for low values of a covariate 
#     so the index is sensitive to ANY presence vs a more steady/linear relationship. E.g. if the curve 
#     is logarithmic any presence of the sound is impactful vs if it is exponential there needs to be a 
#     good amount of the signal to affect the index.
wd = '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/'
indices_df = fread(paste0(wd, 'acoustic_indices_aggregation/averages/site_avg_acoustic_indices.csv'))

# read in abgqi data
abgqi_df = fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI.csv')) %>%
  mutate(ARU = substr(site, 4,5),
         ARUdevice = substr(site, 4,8)) %>%
  select(-contains('var')) %>%
  select(site, wavs, contains('mean'), ARU, -Unidentified_mean) %>%
  rename(Anthropophony = Anthropophony_mean, Biophony = Biophony_mean, Geophony = Geophony_mean,
         Quiet = Quiet_mean, Interference = Interference_mean)

# Clean sites that are known problems
abgqi_df_no_error = abgqi_df %>%
  filter(wavs >= (144 * 2.5)) # 44 sites

# sites with static
# s2lam050_190411, s2lam049_210501 : a lot of static but still has other signal
error_sites = read.csv('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/data/identified_problem_sites_witherrors.csv')
error_sites = paste0(error_sites$SiteID, collapse = '|')
abgqi_df_no_error = abgqi_df_no_error[!grepl(error_sites, abgqi_df_no_error$site),] # 7 sites dropped after filter above

# Number of minutes in model data
mod_df = abgqi_df_no_error

# scale covariates ABGQI
mod_df = mod_df %>%
  mutate(Anthropophony = min_max_norm(Anthropophony),
         Biophony      = min_max_norm(Biophony),
         Geophony      = min_max_norm(Geophony),
         Quiet         = min_max_norm(Quiet),
         Interference  = min_max_norm(Interference))

# Read in model objects
model_objects = readRDS(file = paste0(wd, 'modeling/acoustic_indices/LM_model_objects/gam_model_objects_20June2022.RData'))
model_slopes = readRDS(file = paste0(wd, 'modeling/acoustic_indices/LM_model_objects/gam_model_slopes_20June2022.RData'))



######################################
# Data summaries
length(unique(mod_df$site)) # sites used in analyses
sum(mod_df$wavs)            # number of recordings in analyses

# count number of sites per year
mod_df %>%
  mutate(YY = substr(site, 10, 11)) %>%
  group_by(YY) %>%
  count()

# get number of recordings per year
mod_df %>%
  mutate(YY = substr(site, 10, 11)) %>%
  group_by(YY) %>%
  summarise(sum(wavs))

length(unique(abgqi_df$site)) # total sites pre-analysis
sum(abgqi_df$wavs)            # total wavs pre-analysis

# get mean and sd ABGQ after removing Int
mean_sd <- list(
  mean = ~mean(.x, na.rm = TRUE), 
  sd = ~sd(.x, na.rm = TRUE)
)
mod_df %>%
  select(Anthropophony, Biophony, Geophony, Quiet) %>%
  dplyr::summarise(across(where(is.numeric), mean_sd)) %>%
  round(6) * 100

######################################
# Visualize slope objects
# calculate all slope values
slope_df = do.call("rbind", model_slopes)

# use middle 99% of x values to plot slopes to avoid extreme tail behavior
abgqi_limits = mod_df %>%
  select(Anthropophony, Biophony, Geophony, Quiet, Interference) %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(lower = quantile(value, 0.005),
            upper = quantile(value, 0.995))

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
  group_by(index, var) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) %>%
  mutate(index = factor(index))

# summarize slope trends using CI of mean
# slope_cis = slope_df %>%
#   select(-smooth, -se, -crit, -lower, -upper) %>%
#   group_by(index, var) %>%
#   summarise(lower025 = quantile(derivative, 0.025),
#             upper975 = quantile(derivative, 0.975),
#             median = quantile(derivative, 0.50)) %>%
#   mutate(index = factor(index))

# Facet = soundscape components
# simple 95% CI error plot of summary
avg_slope_error_plot = function(slope_avg_df, temp_index) {
  temp_min = floor(min(slope_avg_df$lower))
  temp_max = ceiling(max(slope_avg_df$upper))
  gg = slope_avg_df %>%
    group_by(var) %>%
    filter(var == temp_index) %>%
    mutate(index = fct_reorder(index, derivative)) %>%
    ggplot(aes(x = index, y = derivative)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ylim(c(temp_min, temp_max)) + # min max in slope df
    ylab("Slope") +
    xlab("Acoustic Index") +
    coord_flip() +
    facet_wrap(~var) +
    theme_bw()
  return(gg)
}

ggAnthro = avg_slope_error_plot(slope_cis, 'Anthropophony')
ggBio = avg_slope_error_plot(slope_cis, 'Biophony')
ggGeo = avg_slope_error_plot(slope_cis, 'Geophony')
ggQuiet = avg_slope_error_plot(slope_cis, 'Quiet')
ggInt = avg_slope_error_plot(slope_cis, 'Interference')
allplots = ggarrange(ggAnthro, ggBio, ggGeo, ggQuiet, ggInt, 
                     ncol = 3, nrow = 2,
                     labels = c("A", "B", "C", "D", "E"))
annotate_figure(allplots, top = text_grob("95% CI for first derivative (slopes) of GAM partial effects.
Note: should be interpreted alongside partial effects plots.", size = 14, face = "bold"),
                bottom = text_grob("Values reflect middle 99% of covariate range based on erroneous tail behavior.",
                                   size = 10))          

# Figure R1
# Facet = GAM models
temp_min = floor(min(slope_cis$lower))
temp_max = ceiling(max(slope_cis$upper))
(gg = slope_cis %>%
  mutate(var = factor(var, levels = c('Interference', 'Quiet','Geophony', 'Biophony',"Anthropophony"))) %>%
  mutate(index = factor(index, levels = c("zcr_mean","ACI","M","H","Hs","NDSI_A","AEI","Ht","rugo","NDSI_B","NDSI",
                             "ADI","sfm","BI","R"))) %>% # ordered by deviance 
  ggplot(aes(x = var, y = derivative)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(temp_min, temp_max)) + # min max in slope df
  ylab("Slope") +
  xlab("Acoustic Index") +
  coord_flip() +
  facet_wrap(~index) +
  theme_bw())

(gg = annotate_figure(gg, top = text_grob("95% CI for first derivative (slopes) of GAM partial effects.
Note: should be interpreted alongside partial effects plots.", size = 14, face = "bold"),
                bottom = text_grob("Values reflect middle 99% of covariate range to minimize erroneous tail behavior.",
                                   size = 10)))

ggsave(gg, filename = paste0(wd, '/modeling/acoustic_indices/GAM_slopes_by_model.png'), 
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
# y-axis is the centered smooth
draw(model_objects$zcr_mean, scales = "fixed")
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

draw(model_objects$M, scales = "fixed")
# minimal biophony effect
# strong log shaped Anthro and exponential Interference
# Geophony is positive until some mixed behavior at sparser mid-high values

draw(model_objects$H, scales = "fixed")

draw(model_objects$Hs, scales = "fixed")

draw(model_objects$NDSI_A, scales = "fixed")

draw(model_objects$AEI, scales = "fixed")

draw(model_objects$Ht, scales = "fixed")

draw(model_objects$rugo, scales = "fixed")

draw(model_objects$NDSI_B, scales = "fixed")

draw(model_objects$NDSI, scales = "fixed")
model_objects$NDSI
# weakly impacted by interference (positive)
# negative Anthro
# Strong positive effect from Bio and lessser so Quiet

draw(model_objects$ADI, scales = "fixed")

draw(model_objects$sfm, scales = "fixed")
model_objects$sfm
summary(indices_df$sfm)

draw(model_objects$BI, scales = "fixed")

draw(model_objects$R, scales = "fixed")
