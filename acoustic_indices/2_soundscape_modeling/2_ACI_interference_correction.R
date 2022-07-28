# 
# 
# # To Do:
# 1) filter for Biophony sites
# 2) filter for biophony and anthro sites
# 3) filter for biophony and interference sites

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
source('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/code/paper1-AcousticIndices/utility_fxs.R')

##########################################
# Read in by-min abgqi data and remove recordings with interference (Brute Force)
results_dir = '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/ABGQI_inference/'

# SITES WITH ANOMALOUS BEHAVIOR 
error_sites = c('s2lam012_190506.csv','s2lam012_190529.csv','s2lam012_190605.csv',
                's2lam018_190602.csv','s2llg001_170606.csv','s2llg002_170418.csv',
                's2llg002_170718.csv','s2llg002_180330.csv','s2llg002_180414.csv',
                's2llg002_180517.csv','s2llg002_180522.csv','s2llg002_190521.csv',
                's2llg002_190530.csv','s2llg004_180330.csv','s2llg006_170627.csv',
                's2llg006_17062.csv')

site_csvs = list.files(paste0(results_dir,"by_site"), full.names = F, pattern = "*.csv")

# Remove any of the above error sites
site_csvs = site_csvs[!site_csvs %in% error_sites]

# function to pull desired stats : will create new col with original Col_name + "mean" or "var"
mean_var <- list(
  mean = ~mean(.x, na.rm = TRUE) 
)

# instantiate empty lists
site_avg_list = list()
wav_list = list()
for(i in 1:length(site_csvs)){
  temp_site = site_csvs[i]
  temp_site_name = strsplit(temp_site, ".csv")[[1]][1]
  print(paste0(i," : ", temp_site_name))
  
  #read in csv
  df = read.csv(paste0(results_dir,"by_site/",temp_site))
  n_wavs = length(unique(df$wav)) # number of recordings
  
  # get by minute presence of interference
  df_int = df %>%
    filter(Interference == 1)
  int_wavs = unique(df_int$wav)
  
  # filter out ANY recording with interference (i.e., >2 sec)
  df_sans_int = df[!df$wav %in% int_wavs,]
  
  ##### SITE AVG #####
  # get column-wise statistics on integer columns only (e.g. binarized)
  temp_avg = df_sans_int %>%
    dplyr::summarise(across(where(is.integer), mean_var))

  # list of wavs 
  wav_list[[temp_site_name]] = unique(tools::file_path_sans_ext(df_sans_int$wav))
  
  # store site name and stats
  site_avg_list[[i]] = as.data.frame(c("site" = temp_site_name, 
                                       "wavs" = n_wavs, 
                                       "wavs_sans_int" = length(int_wavs),
                                       temp_avg))
}

# concat all site avgs
all_site_avg_df = do.call("rbind", site_avg_list)
wavs_sans_int = data.frame("no_int" = unlist(wav_list), row.names = NULL)

# save csv
write.csv(all_site_avg_df, file = paste0(results_dir, "/averages/site_avg_ABGQI_sans_interference.csv"), row.names = FALSE)
write.csv(wavs_sans_int, file = paste0(results_dir, "/averages/wavs_sans_interference.csv"), row.names = FALSE)

###########################################
####### filter acoustic indices  ##########
results_dir = '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/acoustic_indices_aggregation/'

# SITES WITH ANOMALOUS BEHAVIOR 
error_sites = c('s2lam012_190506.csv','s2lam012_190529.csv','s2lam012_190605.csv',
                's2lam018_190602.csv','s2llg001_170606.csv','s2llg002_170418.csv',
                's2llg002_170718.csv','s2llg002_180330.csv','s2llg002_180414.csv',
                's2llg002_180517.csv','s2llg002_180522.csv','s2llg002_190521.csv',
                's2llg002_190530.csv','s2llg004_180330.csv','s2llg006_170627.csv',
                's2llg006_17062.csv')


# 1-min site acoustic indices
site_csvs = list.files(paste0(results_dir,"by_site"), full.names = F, pattern = "*.csv")

# Remove any of the above error sites
site_csvs = site_csvs[!site_csvs %in% error_sites]

# instantiate empty lists
site_avg_list = list()
for(i in 1:length(site_csvs)){
  temp_site = site_csvs[i]
  temp_site_name = strsplit(temp_site, ".csv")[[1]][1]
  print(paste0(i," : ", temp_site_name))
  
  #read in csv
  df = read.csv(paste0(results_dir,"by_site/",temp_site)) 

  # filter based on wavs with no interference
  temp_wav_vec = wav_list[[temp_site_name]]
  df_sans_int = df[df$wav %in% temp_wav_vec,]
  
  ##### SITE AVG #####
  # get column-wise statistics on integer columns only (e.g. binarized)
  temp_avg = df_sans_int %>%
    select(-c(YYYY, MM, DD, hh, mm, zcr_max, zcr_min)) %>%
    dplyr::summarise(across(where(is.numeric), mean))

  # store site name and stats
  site_avg_list[[i]] = as.data.frame(c("site" = temp_site_name, temp_avg))
  
}

# concat all site avgs
all_site_avg_df = do.call("rbind", site_avg_list)

# save csv
write.csv(all_site_avg_df, file = paste0(results_dir, "/averages/site_avg_acoustic_indices_sans_interference.csv"), row.names = FALSE)


###########################################
####### read in data and data prep ########
wd = '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/'

# read in abgqi data
abgqi_df = fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI_sans_interference.csv')) %>%
  mutate(ARU = substr(site, 4,5),
         ARUdevice = substr(site, 4,8)) %>%
  select(site, wavs, wavs_sans_int, contains('mean'), ARU, -Interference_mean) %>%
  rename(Anthropophony = Anthropophony_mean, Biophony = Biophony_mean, Geophony = Geophony_mean,
         Quiet = Quiet_mean)

# Clean sites that are known problems
abgqi_df_no_error = abgqi_df %>%
  filter(wavs >= (144 * 2.5)) # 44 sites

# sites with static
# s2lam050_190411, s2lam049_210501 : a lot of static but still has other signal
error_sites = read.csv('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/data/identified_problem_sites_witherrors.csv')
error_sites = paste0(error_sites$SiteID, collapse = '|')
abgqi_df_no_error = abgqi_df_no_error[!grepl(error_sites, abgqi_df_no_error$site),] # 7 sites dropped after filter above

indices_df = fread(paste0(wd, 'acoustic_indices_aggregation/averages/site_avg_acoustic_indices_sans_interference.csv'))

# Number of minutes in model data
mod_df = abgqi_df_no_error

# scale covariates ABGQI
mod_df = mod_df %>%
  filter(complete.cases(.)) %>%
  mutate(Anthropophony = min_max_norm(Anthropophony),
         Biophony      = min_max_norm(Biophony),
         Geophony      = min_max_norm(Geophony),
         Quiet         = min_max_norm(Quiet))

sum(mod_df$wavs_sans_int)

# get mean and sd ABGQ after removing Int
mean_sd <- list(
  mean = ~mean(.x, na.rm = TRUE), 
  sd = ~sd(.x, na.rm = TRUE)
)
mod_df %>%
  select(Anthropophony, Biophony, Geophony, Quiet) %>%
  dplyr::summarise(across(where(is.numeric), mean_sd)) %>%
  round(6)

###########################################
################# ACI #####################
i = "ACI"
temp_y = data.frame(site = indices_df$site, ACI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>%
  mutate(ARU = factor(ARU),
         logwavs = log(wavs))

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
# bring ACI down to 0
temp_df$ACI_0 = temp_df$ACI - min(temp_df$ACI) + 1e-7

# both ARUs have similar shape and wiggliness -- if Randeffect is included use GS model type
ggplot(temp_df, aes(x = ACI_0, fill = ARU)) +
  geom_histogram()

mod1 = gam(ACI_0 ~
             ARU +
             logwavs +
             s(Anthropophony, k = 5) +
             s(Biophony, k = 5) +
             s(Geophony, k = 5) +
             s(Quiet, k = 5),
           data = temp_df,
           family = Gamma(link = 'log'),
           method = 'ML')
summary(mod1)
par(mfrow = c(2, 2))
gam.check(mod1) # heavy tail in QQ - primarily de to one observation

resids = data.frame(mod1$residuals) # most extreme obs is the only ACI = 0
temp_df_remod = temp_df[-67,]

# remodel without extreme outlier
remod1 = gam(ACI_0 ~
               ARU +
               logwavs +
               s(Anthropophony, k = 5) +
               s(Biophony, k = 5) +
               s(Geophony, k = 5) +
               s(Quiet, k = 5),
             data = temp_df_remod,
             family = Gamma(link = 'log'),
             method = 'ML')
summary(remod1)
par(mfrow = c(2, 2))
gam.check(remod1, rep = 500) # good diagnostics now
AIC(mod1, remod1)

# remove logwavs
mod2 = gam(ACI_0 ~
             ARU +
             s(Anthropophony, k = 5) +
             s(Biophony, k = 5) +
             s(Geophony, k = 5) +
             s(Quiet, k = 5),
           data = temp_df_remod,
           family = Gamma(link = 'log'),
           method = 'ML')
summary(mod2)
gam.check(mod2, rep = 500)
AIC(remod1, mod2)

# remove Quiet
mod3 = gam(ACI_0 ~
             ARU +
             s(Anthropophony, k = 5) +
             s(Biophony, k = 5) +
             s(Geophony, k = 5),
           data = temp_df_remod,
           family = Gamma(link = 'log'),
           method = 'ML')
summary(mod3)
AIC(mod2, mod3) # Quiet suggested to stay in
gam.check(mod3, rep = 500)


# visualize partial effects
draw(mod2, scales = "fixed", residuals = TRUE)
appraise(mod2, method = "simulate")
summary(mod2)


aci_with_int = readRDS(file = paste0(wd, 'modeling/acoustic_indices/LM_model_objects/gam_model_objects_20June2022.RData'))$ACI
draw(aci_with_int, scales = "fixed")
summary(aci_with_int)
##### Generate Slope figure 
# summarized slopes
ci = slope_summary(mod2)
model_slopes = cbind(index = rep(i, times = nrow(ci)), ci)

# visualize predictions vs observed
pred_plot(data_df = temp_df, model_fit = mod2, index_name = "ACI_0")

# use middle 99% of x values to plot slopes to avoid extreme tail behavior
abgqi_limits = mod_df %>%
  select(Anthropophony, Biophony, Geophony, Quiet) %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(lower = quantile(value, 0.005),
            upper = quantile(value, 0.995))

# filter slope df based on percentile values
slope_df_filtered = model_slopes %>%
  filter(case_when(var == 'Anthropophony' ~ data >= abgqi_limits$lower[1] & data <= abgqi_limits$upper[1],
                   var == 'Biophony'      ~ data >= abgqi_limits$lower[2] & data <= abgqi_limits$upper[2],
                   var == 'Geophony'      ~ data >= abgqi_limits$lower[3] & data <= abgqi_limits$upper[3],
                   var == 'Quiet'         ~ data >= abgqi_limits$lower[4] & data <= abgqi_limits$upper[4]))


# summarize slope trends
slope_cis = slope_df_filtered %>%
  select(-smooth) %>%
  group_by(index, var) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) %>%
  mutate(index = factor(index))

######
# Interference model
model_slopes_int = readRDS(file = paste0(wd, 'modeling/acoustic_indices/LM_model_objects/gam_model_slopes_20June2022.RData'))$ACI
abgqi_df_int = fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI.csv')) %>%
  mutate(ARU = substr(site, 4,5),
         ARUdevice = substr(site, 4,8)) %>%
  select(-contains('var')) %>%
  select(site, wavs, contains('mean'), ARU, -Unidentified_mean) %>%
  rename(Anthropophony = Anthropophony_mean, Biophony = Biophony_mean, Geophony = Geophony_mean,
         Quiet = Quiet_mean, Interference = Interference_mean)

# Clean sites that are known problems
abgqi_df_no_error_int = abgqi_df_int %>%
  filter(wavs >= (144 * 2.5)) # 44 sites
mod_df_int = abgqi_df_no_error_int[!grepl(error_sites, abgqi_df_no_error_int$site),]

abgqi_limits_int = mod_df_int %>%
  select(Anthropophony, Biophony, Geophony, Quiet, Interference) %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(lower = quantile(value, 0.005),
            upper = quantile(value, 0.995))
slope_df_filtered_int = model_slopes_int %>%
  filter(case_when(var == 'Anthropophony' ~ data >= abgqi_limits_int$lower[1] & data <= abgqi_limits_int$upper[1],
                   var == 'Biophony'      ~ data >= abgqi_limits_int$lower[2] & data <= abgqi_limits_int$upper[2],
                   var == 'Geophony'      ~ data >= abgqi_limits_int$lower[3] & data <= abgqi_limits_int$upper[3],
                   var == 'Interference'  ~ data >= abgqi_limits_int$lower[4] & data <= abgqi_limits_int$upper[4],
                   var == 'Quiet'         ~ data >= abgqi_limits_int$lower[5] & data <= abgqi_limits_int$upper[5]))
slope_cis_int = slope_df_filtered_int %>%
  select(-smooth) %>%
  group_by(index, var) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) %>%
  mutate(index = factor(index))

slope_cis$Int = "no Interference"
slope_cis_int$Int = "Interference"
slope_cis = rbind(slope_cis, slope_cis_int)

temp_min = floor(min(slope_cis$lower))
temp_max = ceiling(max(slope_cis_int$upper))
# Plot int and no-int slope CIs
(gg = slope_cis %>%
    mutate(var = factor(var, levels = c('Interference','Quiet','Geophony', 'Biophony','Anthropophony'))) %>%
    ggplot(aes(x = var, y = derivative)) +
    geom_errorbar(aes(ymin = lower, ymax = upper, linetype = factor(Int)), 
                  width = 0.3, size = 0.7, position = position_dodge()) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ylim(c(temp_min, temp_max)) + # min max in slope df
    ylab("Slope") +
    xlab("Acoustic Index") +
    coord_flip() +
    facet_wrap(~index) +
    guides(linetype = guide_legend(title = "GAM variant")) +
    theme_bw())

(gg = annotate_figure(gg, top = text_grob("95% CI for first derivative (slopes) of GAM partial effects.
Note: should be interpreted alongside partial effects plots.", size = 14, face = "bold"),
                      bottom = text_grob("Values reflect middle 99% of covariate range to minimize erroneous tail behavior.",
                                         size = 10)))
