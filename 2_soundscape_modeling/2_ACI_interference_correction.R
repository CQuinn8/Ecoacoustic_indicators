library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(mgcv)
library(gratia)

# Lib that may help install fonts if errors occur
# library(remotes)
# remotes::install_version("Rttf2pt1", version = "1.3.8")

library(extrafont)
# font_import() # run once if needed
loadfonts(device = "win")
fonts()

# custom additional functions
source('utility_fxs.R')

##########################################
# Read in by-min abgqi data and remove recordings with interference (Brute Force)
results_dir <- '/paper-AcousticIndices/results/ABGQI_inference/'

# SITES WITH ANOMALOUS BEHAVIOR 
error_sites <- c('s2lam012_190506.csv','s2lam012_190529.csv','s2lam012_190605.csv',
                's2lam018_190602.csv','s2llg001_170606.csv','s2llg002_170418.csv',
                's2llg002_170718.csv','s2llg002_180330.csv','s2llg002_180414.csv',
                's2llg002_180517.csv','s2llg002_180522.csv','s2llg002_190521.csv',
                's2llg002_190530.csv','s2llg004_180330.csv','s2llg006_170627.csv',
                's2llg006_17062.csv')

# List all sites  csvs
site_csvs <- list.files(paste0(results_dir, "by_site"), full.names = F, pattern = "*.csv")

# Remove any of the above error sites
site_csvs <- site_csvs[!site_csvs %in% error_sites]

# function to pull desired stats : will create new col with original Col_name + "mean" or "var" etc
mean_fx <- list(
  mean = ~mean(.x, na.rm = TRUE) 
)

# save output csv of avg values without interference
out_file <- paste0(results_dir, "/averages/site_avg_ABGQI_sans_interference.csv")
if(file.exists(out_file)){
  abgqi_df <- read.csv(out_file)
} else {
  # instantiate empty lists
  site_avg_list <- list()
  wav_list <- list()
  
  # Iterate over every site's by-min ABGQI pres/abs 
  for(i in 1:length(site_csvs)){
    temp_site <- site_csvs[i]
    temp_site_name <- strsplit(temp_site, ".csv")[[1]][1]
    print(paste0(i," : ", temp_site_name))
    
    #read in csv
    df <- read.csv(paste0(results_dir,"by_site/",temp_site))
    
    # number of recordings
    n_wavs <- length(unique(df$wav))
    
    # get any recordings with interference
    df_int <- df %>%
      filter(Interference == 1)
    int_wavs <- unique(df_int$wav)
    
    # filter out ANY recording with interference (i.e., >2 sec)
    df_sans_int <- df[!df$wav %in% int_wavs,]
    
    ##### SITE AVG #####
    # get column-wise statistics on integer columns only
    temp_avg <- df_sans_int %>%
      dplyr::summarise(across(where(is.integer), mean_fx))
  
    # list of wavs 
    wav_list[[temp_site_name]] <- unique(tools::file_path_sans_ext(df_sans_int$wav))
    
    # store site name and stats
    site_avg_list[[i]] <- as.data.frame(c("site" = temp_site_name, 
                                         "wavs" = n_wavs, 
                                         "wavs_sans_int" = length(int_wavs),
                                         temp_avg))
  }
  
  # concat all site avgs
  all_site_avg_df <- do.call("rbind", site_avg_list)
  wavs_sans_int <- data.frame("no_int" = unlist(wav_list), row.names = NULL)
  
  # save csv
  write.csv(all_site_avg_df, file = out_file, row.names = FALSE)
  write.csv(wavs_sans_int, file = paste0(results_dir, "/averages/wavs_sans_interference.csv"), row.names = FALSE)
}
rm(results_dir)

###########################################
####### filter acoustic indices  ##########
# This section uses the above dataframe where Interference containing recordings have been rm
results_dir <- '/paper-AcousticIndices/results/acoustic_indices_aggregation/'

# 1-min site acoustic indices
site_csvs <- list.files(paste0(results_dir, "by_site"), full.names = F, pattern = "*.csv")

# Remove any of the above error sites
site_csvs <- site_csvs[!site_csvs %in% error_sites]

out_file <- paste0(results_dir, "/averages/site_avg_acoustic_indices_sans_interference.csv")
if(file.exists(out_file)){
  indices_df <- read.csv(out_file)
} else {
  site_avg_list <- list()
  
  # Read over every acoustic index site csv
  for(i in 1:length(site_csvs)){
    temp_site <- site_csvs[i]
    temp_site_name <- strsplit(temp_site, ".csv")[[1]][1]
    print(paste0(i," : ", temp_site_name))
    
    #read in csv
    df <- read.csv(paste0(results_dir,"by_site/",temp_site)) 
  
    # filter based on wavs with no interference
    temp_wav_vec <- wav_list[[temp_site_name]]
    df_sans_int <- df[df$wav %in% temp_wav_vec,]
    
    ##### SITE AVG #####
    # get column-wise statistics on integer columns only
    temp_avg <- df_sans_int %>%
      select(-c(YYYY, MM, DD, hh, mm, zcr_max, zcr_min)) %>%
      dplyr::summarise(across(where(is.numeric), mean))
  
    # store site name and stats
    site_avg_list[[i]] <- as.data.frame(c("site" = temp_site_name, temp_avg))
    
  }
  
  # concat all site avgs
  all_site_avg_df <- do.call("rbind", site_avg_list)
  
  # save csv
  write.csv(all_site_avg_df, file = out_file, row.names = FALSE)
}


###########################################
####### read in data and data prep ########
# Code from here on is similar to 0_by_site_GAM.R
wd <- '/paper-AcousticIndices/results/'

# read in abgqi data
abgqi_df <- fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI_sans_interference.csv')) %>%
  mutate(ARU = substr(site, 4,5),
         ARUdevice = substr(site, 4,8)) %>%
  select(site, wavs, wavs_sans_int, contains('mean'), ARU, -Interference_mean) %>%
  rename(Anthropophony = Anthropophony_mean, Biophony = Biophony_mean, Geophony = Geophony_mean,
         Quiet = Quiet_mean)

# Clean sites that are known problems
abgqi_df_no_error <- abgqi_df %>%
  filter(wavs >= (144 * 2.5)) # 44 sites

# sites with static
error_sites <- read.csv('/paper-AcousticIndices/data/identified_problem_sites_witherrors.csv')
error_sites <- paste0(error_sites$SiteID, collapse = '|')
abgqi_df_no_error <- abgqi_df_no_error[!grepl(error_sites, abgqi_df_no_error$site),] # 7 sites dropped after filter above

# Read in the updated, interference filtered acoustic index csv
indices_df <- fread(paste0(wd, 'acoustic_indices_aggregation/averages/site_avg_acoustic_indices_sans_interference.csv'))

# scale covariates ABGQ
mod_df <- abgqi_df_no_error %>%
  filter(complete.cases(.)) %>%
  mutate(Anthropophony = min_max_norm(Anthropophony),
         Biophony      = min_max_norm(Biophony),
         Geophony      = min_max_norm(Geophony),
         Quiet         = min_max_norm(Quiet))

# count of recordings without interference (n = 300,157)
sum(mod_df$wavs_sans_int)

# Summarize, non-interference ABGQ
mod_df %>%
  select(Anthropophony, Biophony, Geophony, Quiet) %>%
  dplyr::summarise(across(where(is.numeric), mean_sd)) %>%
  round(6)


###########################################
######## ACI - sans Interference ##########
i <- "ACI"
temp_y <- data.frame(site = indices_df$site, ACI = indices_df[[i]])
temp_df <- merge(x = temp_y, y = mod_df) %>%
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
temp_df$ACI_0 <- temp_df$ACI - min(temp_df$ACI) + 1e-7
ggplot(temp_df, aes(x = ACI_0, fill = ARU)) +
  geom_histogram()

mod1 <- gam(ACI_0 ~
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

resids <- data.frame(mod1$residuals) # most extreme obs is the only ACI = 0
temp_df_remod <- temp_df[-67,]

# remodel without extreme outlier
remod1 <- gam(ACI_0 ~
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
mod2 <- gam(ACI_0 ~
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
mod3 <- gam(ACI_0 ~
             ARU +
             s(Anthropophony, k = 5) +
             s(Biophony, k = 5) +
             s(Geophony, k = 5),
           data = temp_df_remod,
           family = Gamma(link = 'log'),
           method = 'ML')
summary(mod3)
mod2$aic - mod3$aic # Quiet suggested to stay in
gam.check(mod3, rep = 500)


# visualize partial effects for final model (MOD2)
appraise(mod2, method = "simulate")
summary(mod2)

# combined visual
(gg <- draw(mod2, scales = "fixed", residuals = TRUE, ncol = 4) &
    theme_bw())

# individual plotting of covariated in order for proper formatting
sm <- smooth_estimates(mod2) %>%
  add_confint()
temp_df_remod <- temp_df_remod %>%
  add_partial_residuals(mod2)

# Repeated GG plot statement for each covariate PDP
(gg1 <- sm %>%
  filter(smooth == "s(Anthropophony)") %>%
  ggplot() +
    geom_rug(aes(x = Anthropophony), 
             data = temp_df_remod,
             sides = "b", length = grid::unit(0.02, "npc")) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Anthropophony), alpha = 0.2) +
    geom_point(aes(x = Anthropophony, y = `s(Anthropophony)`),
               data = temp_df_remod, cex = 1.5, colour = "steelblue3", alpha = 0.1) +
    geom_line(aes(x = Anthropophony, y = est), lwd = 1.1, alpha = 0.8) +
    ylim(-4,4) +
    theme_bw() +
    labs(y = "Partial effect", title = "Anthropophony"))
(gg2 <- sm %>%
    filter(smooth == "s(Biophony)") %>%
    ggplot() +
    geom_rug(aes(x = Biophony), 
             data = temp_df_remod,
             sides = "b", length = grid::unit(0.02, "npc")) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Biophony), alpha = 0.2) +
    geom_point(aes(x = Biophony, y = `s(Biophony)`),
               data = temp_df_remod, cex = 1.5, colour = "steelblue3", alpha = 0.1) +
    geom_line(aes(x = Biophony, y = est), lwd = 1.1, alpha = 0.8) +
    ylim(-4,4) +
    theme_bw() +
    labs(y = "", title = "Biophony"))
(gg3 <- sm %>%
    filter(smooth == "s(Geophony)") %>%
    ggplot() +
    geom_rug(aes(x = Geophony), 
             data = temp_df_remod,
             sides = "b", length = grid::unit(0.02, "npc")) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Geophony), alpha = 0.2) +
    geom_point(aes(x = Geophony, y = `s(Geophony)`),
               data = temp_df_remod, cex = 1.5, colour = "steelblue3", alpha = 0.1) +
    geom_line(aes(x = Geophony, y = est), lwd = 1.1, alpha = 0.8) +
    ylim(-4,4) +
    theme_bw() +
    labs(y = "", title = "Geophony"))
(gg4 <- sm %>%
    filter(smooth == "s(Quiet)") %>%
    ggplot() +
    geom_rug(aes(x = Quiet), 
             data = temp_df_remod,
             sides = "b", length = grid::unit(0.02, "npc")) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Quiet), alpha = 0.2) +
    geom_point(aes(x = Quiet, y = `s(Quiet)`),
               data = temp_df_remod, cex = 1.5, colour = "steelblue3", alpha = 0.1) +
    geom_line(aes(x = Quiet, y = est), lwd = 1.1, alpha = 0.8) +
    ylim(-4,4) +
    theme_bw() +
    labs(y = "", title = "Quiet"))

# Combine PDPs into one plot obj
gg1 + gg2 + gg3 + gg4 + patchwork::plot_layout(ncol = 4)

# Save output (Figure 4)
ggsave(filename = 'figure_4.png', 
       plot = gg1 + gg2 + gg3 + gg4 + patchwork::plot_layout(ncol = 4), 
       device = 'png',
       path = 'figures/',
       width = 8, height = 4, dpi = 500)



###########################################
# Compare with and without Interference ACI
# Read in ACI GAM object with Interference included
aci_with_int <- readRDS(file = '/models/gam_model_objectsmean_fx.RData')$ACI
draw(aci_with_int, scales = "fixed")
summary(aci_with_int)

# Generate slope summaries for sans-Interference
ci <- slope_summary(mod2)
model_slopes <- cbind(index = rep(i, times = nrow(ci)), ci)

# visualize predictions vs observed
pred_plot(data_df = temp_df, model_fit = mod2, index_name = "ACI_0")

# use middle 99% of x values to plot slopes to avoid extreme tail behavior
abgqi_limits <- mod_df %>%
  select(Anthropophony, Biophony, Geophony, Quiet) %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(lower = quantile(value, 0.005),
            upper = quantile(value, 0.995))

# filter slope df based on values corresponding to 99% percentile values
slope_df_filtered <- model_slopes %>%
  filter(case_when(var == 'Anthropophony' ~ data >= abgqi_limits$lower[1] & data <= abgqi_limits$upper[1],
                   var == 'Biophony'      ~ data >= abgqi_limits$lower[2] & data <= abgqi_limits$upper[2],
                   var == 'Geophony'      ~ data >= abgqi_limits$lower[3] & data <= abgqi_limits$upper[3],
                   var == 'Quiet'         ~ data >= abgqi_limits$lower[4] & data <= abgqi_limits$upper[4]))


# PROCESS: Interference model and data
model_slopes_int <- readRDS(file = paste0(wd, '/models/gam_model_slopesmean_fx.RData'))$ACI
abgqi_df_int <- fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI.csv')) %>%
  mutate(ARU = substr(site, 4,5),
         ARUdevice = substr(site, 4,8)) %>%
  select(-contains('var')) %>%
  select(site, wavs, contains('mean'), ARU, -Unidentified_mean) %>%
  rename(Anthropophony = Anthropophony_mean, Biophony = Biophony_mean, Geophony = Geophony_mean,
         Quiet = Quiet_mean, Interference = Interference_mean)

# Clean sites that are known problems
abgqi_df_no_error_int <- abgqi_df_int %>%
  filter(wavs >= (144 * 2.5)) 
mod_df_int <- abgqi_df_no_error_int[!grepl(error_sites, abgqi_df_no_error_int$site),]

abgqi_limits_int <- mod_df_int %>%
  select(Anthropophony, Biophony, Geophony, Quiet, Interference) %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(lower = quantile(value, 0.005),
            upper = quantile(value, 0.995))

slope_df_filtered_int <- model_slopes_int %>%
  filter(case_when(var == 'Anthropophony' ~ data >= abgqi_limits_int$lower[1] & data <= abgqi_limits_int$upper[1],
                   var == 'Biophony'      ~ data >= abgqi_limits_int$lower[2] & data <= abgqi_limits_int$upper[2],
                   var == 'Geophony'      ~ data >= abgqi_limits_int$lower[3] & data <= abgqi_limits_int$upper[3],
                   var == 'Interference'  ~ data >= abgqi_limits_int$lower[4] & data <= abgqi_limits_int$upper[4],
                   var == 'Quiet'         ~ data >= abgqi_limits_int$lower[5] & data <= abgqi_limits_int$upper[5]))


# COMPARE MODELS
# Plot int and no-int slope boxplots
slope_df_filtered$Interference <- factor("No Interference")
slope_df_filtered_int$Interference <- factor("Interference")
slope_df_filtered_combined <- rbind(slope_df_filtered, slope_df_filtered_int)
temp_min <- floor(min(slope_df_filtered_combined$lower))
temp_max <- ceiling(max(slope_df_filtered_combined$upper))

(gg <- slope_df_filtered_combined %>%
    mutate(var = factor(var, levels = c('Interference', 'Quiet', 'Geophony', 'Biophony', 'Anthropophony'))) %>%
    ggplot(aes(x = factor(var), y = derivative, fill = Interference)) +
      geom_boxplot(width = 0.5, 
                   outlier.shape = NA, 
                   position = position_dodge(preserve = "single"),
                   notch = TRUE) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ylim(c(temp_min, temp_max)) + # min max in slope df
      ylab("Slope") +
      xlab("") +
      scale_fill_manual(values = c('#B8B8B8', '#FFFFFF'), name = "") +
      coord_flip() +
      theme_bw() +
      theme(legend.position = 'bottom',
            text = element_text(size = 16, family = 'Calibri'),
            axis.text.y = element_text(angle = 45, vjust = -0.5)))


# Save Fig3
ggsave(filename = 'figure_3.png', 
       plot = gg, 
       device = 'png',
       path = '/figures/',
       width = 4.5, height = 6, dpi = 500)
