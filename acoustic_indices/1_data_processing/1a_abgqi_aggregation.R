# Purpose : Use site csvs of all melspec ABGQI preds to generate by-in, by-hour, and by-site averages

ll = '/projects/tropics/users/cquinn/R_400/'
.libPaths(c( .libPaths(), ll))
library(data.table)
library(dplyr)

results_dir = '/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/ABGQI_inference/'

# SITES WITH ANOMALOUS BEHAVIOR 
error_sites = c('s2lam012_190506.csv','s2lam012_190529.csv','s2lam012_190605.csv',
                's2lam018_190602.csv','s2llg001_170606.csv','s2llg002_170418.csv',
                's2llg002_170718.csv','s2llg002_180330.csv','s2llg002_180414.csv',
                's2llg002_180517.csv','s2llg002_180522.csv','s2llg002_190521.csv',
                's2llg002_190530.csv','s2llg004_180330.csv','s2llg006_170627.csv',
                's2llg006_17062.csv')



#######################################
# 1min by-site Average
#######################################
site_csvs = list.files(paste0(results_dir,"by_site"), full.names = F, pattern = "*.csv")

# Remove any of the above error sites
site_csvs = site_csvs[!site_csvs %in% error_sites]

# function to pull desired stats : will create new col with original Col_name + "mean" or "var"
mean_var <- list(
  mean = ~mean(.x, na.rm = TRUE), 
  var = ~var(.x, na.rm = TRUE)
)

# instantiate empty lists
site_avg_list = list()
by_hour_list = list()
by_hour_all_list = list()
for(i in 1:length(site_csvs)){
  temp_site = site_csvs[i]
  temp_site_name = strsplit(temp_site, ".csv")[[1]][1]
  print(paste0(i," : ", temp_site_name))
  
  #read in csv
  df = read.csv(paste0(results_dir,"by_site/",temp_site))
  n_wavs = length(unique(df$wav)) # number of recordings
    
  # force a sixth class to positive if no other labels are present (akin to unidentifiable)
  df = df %>%
    mutate(Unidentified = ifelse(df$"Anthropophony" == 0 & 
                               df$"Biophony" == 0 & 
                               df$"Geophony" == 0 & 
                               df$"Interference" == 0 &
                               df$"Quiet" == 0,
                             1,0))
  df$Unidentified = as.integer(df$Unidentified)

  ###### SITE BY HOUR (all) ######
  # create DD column based on wav name
  df$DD = substr(df$wav, 25, 26)
  df$MM = substr(df$wav, 22, 23)
  df$HH = substr(df$wav, 28, 29)
  df_by_hour_all = df %>%
    select(-c(wav, mfcc))  
  
  # get column-wise statistics on integer columns only (e.g. binarized)
  temp_hr_all = df_by_hour_all %>%
    group_by(MM, DD, HH) %>%
    dplyr::summarise(across(where(is.numeric), mean))
  
  # count n wavs per hour
  n_all = df_by_hour_all %>%
    group_by(MM,DD,HH) %>%
    dplyr::summarise(n = n() / 30)
  
  # join count and stats 
  temp_hr_avg_all = merge(x = temp_hr_all, y = n_all, by = c('MM','DD','HH'))
    
  
  # store site name and stats
  by_hour_all_list[[i]] = as.data.frame(c("site" = temp_site_name, temp_hr_avg_all))
  
  
  # ###### SITE BY HOUR (24) ######
  # create HH column based on wav name
  df_by_hour = df
  HH_temp = strsplit(as.character(df_by_hour$wav), "_")
  HH_temp = lapply(HH_temp, function(x) x[length(x)]) # pull out final chars ('HH_mm.csv')
  HH_temp = as.integer(lapply(HH_temp, function(x) strsplit(x, "-")[[1]][1])) # pull out HH
  df_by_hour$HH = HH_temp

  # get column-wise statistics on integer columns only (e.g. binarized)
  temp_hr = df_by_hour %>%
    group_by(HH) %>%
    dplyr::summarise(across(where(is.integer), mean_var))

  # count n wavs per hour
  n = df_by_hour %>%
    group_by(HH) %>%
    dplyr::summarise(n = n()/30)

  # join count and stats
  temp_hr_avg = merge(x = temp_hr, y = n, by = 'HH')

  # store site name and stats
  by_hour_list[[i]] = as.data.frame(c("site" = temp_site_name, temp_hr_avg))

  ##### SITE AVG #####
  # get column-wise statistics on integer columns only (e.g. binarized)
  temp_avg = df %>%
    dplyr::summarise(across(where(is.integer), mean_var))
  temp_count = df %>%
    dplyr::summarise(across(where(is.integer), sum))

  # store site name and stats
  site_avg_list[[i]] = as.data.frame(c("site" = temp_site_name, "wavs" = n_wavs, "preds" = 30*n_wavs,
                                  temp_avg, temp_count))
}

# concat all site avgs
all_site_avg_df = do.call("rbind", site_avg_list)
all_site_by_hour_df = do.call("rbind", by_hour_list)
all_site_by_hour_all_df = do.call("rbind", by_hour_all_list)

# save csv
write.csv(all_site_avg_df, file = paste0(results_dir, "/averages/site_avg_ABGQI.csv"), row.names = FALSE)
write.csv(all_site_by_hour_df, file = paste0(results_dir, "/averages/site_by_hour_ABGQI.csv"), row.names = FALSE)
write.csv(all_site_by_hour_all_df, file = paste0(results_dir, "/averages/site_by_hour_all_ABGQI.csv"), row.names = FALSE)
