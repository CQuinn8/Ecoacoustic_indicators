# Purpose : Use site csvs of all acoustic index calculations to generate by-in, by-hour, and by-site data

ll = '/projects/tropics/users/cquinn/R_400/'
.libPaths(c( .libPaths(), ll))
library(data.table)
library(dplyr)

results_dir = '/projects/tropics/users/cquinn/s2l/paper-AcousticIndices/results/acoustic_indices_aggregation/'

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

# dawn chorus (4 am to 12 pm)
dawn_chorus = seq(from = 4, to = 11, by = 1)

# instantiate empty lists
site_avg_list = list()
by_hour_list = list()
by_hour_all_list = list()
dawn_chorus_list = list()
options(dplyr.summarise.inform = FALSE)
for(i in 1:length(site_csvs)){
  temp_site = site_csvs[i]
  temp_site_name = strsplit(temp_site, ".csv")[[1]][1]
  print(paste0(i," : ", temp_site_name))
  
  #read in csv
  df = read.csv(paste0(results_dir,"by_site/",temp_site)) 
  n_wavs = nrow(df) # number of recordings
  
  ###### SITE BY HOUR (all recorded day and hour) ######
  # create DD column based on wav name
  df$DD = substr(df$wav, 25, 26)
  df$MM = substr(df$wav, 22, 23)

  df_by_hour_all = df %>%
    select(-c(wav, mm)) %>%
    rename(HH = hh) 
  
  # get column-wise statistics on integer columns only (e.g. binarized)
  temp_hr_all = df_by_hour_all %>%
    group_by(MM, DD, HH) %>%
    dplyr::summarise(across(where(is.numeric), mean))
  
  # count n wavs per hour
  n_all = df_by_hour_all %>%
    group_by(MM,DD,HH) %>%
    dplyr::summarise(n = n())
  
  # join count and stats 
  temp_hr_avg_all = merge(x = temp_hr_all, y = n_all, by = c('MM','DD','HH'))
  temp_hr_avg_all$YYYY = paste0('20',substr(temp_site_name, 10, 11))
  
  # store site name and stats
  by_hour_all_list[[i]] = as.data.frame(c("site" = temp_site_name, temp_hr_avg_all))
  

  ###### SITE BY HOUR (24) ######
  # create HH column based on wav name
  df_by_hour = df %>%
    select(-c(file, MM, DD, mm, DDhh, hhmm)) %>%
    rename(HH = hh)

  # get column-wise statistics on integer columns only (e.g. binarized)
  temp_hr = df_by_hour %>%
    group_by(HH) %>%
    dplyr::summarise(across(where(is.numeric), mean))

  # count n wavs per hour
  n = df_by_hour %>%
    group_by(HH) %>%
    dplyr::summarise(n = n())

  # join count and stats
  temp_hr_avg = merge(x = temp_hr, y = n, by = 'HH')

  # store site name and stats
  by_hour_list[[i]] = as.data.frame(c("site" = temp_site_name, temp_hr_avg))

  ##### SITE AVG #####
  # get column-wise statistics on integer columns only (e.g. binarized)
  temp_avg = df %>%
    select(-c(file, MM, DD, mm, DDhh, hhmm, hh)) %>%
    dplyr::summarise(across(where(is.numeric), mean))

  # store site name and stats
  site_avg_list[[i]] = as.data.frame(c("site" = temp_site_name, temp_avg))
  
  ###### SITE DAWN CHORUS ######
  if(any(unique(df$hh) %in% dawn_chorus)){
    temp_dawn_avg = df %>% 
      select(-wav, -YYYY, -MM, -DD, -mm) %>%
      filter(hh %in% dawn_chorus) %>%
      summarise(across(where(is.numeric), mean)) %>%
      select(-hh)
    
    # count n wavs per hour
    n = df %>%
      filter(hh %in% dawn_chorus) %>%
      summarise('n_wavs' = n_distinct(wav))
    
    # store site name and stats
    dawn_chorus_list[[i]] = as.data.frame(c("site" = temp_site_name, temp_dawn_avg, 'wavs' = n$n_wavs))
    
  } else {
    print('Site contains no dawn chorus samples')
  }
  
}

# concat all site avgs
all_site_avg_df = do.call("rbind", site_avg_list)
all_site_by_hour_df = do.call("rbind", by_hour_list)
all_site_by_hour_all_df = do.call("rbind", by_hour_all_list)
all_site_dawn_chorus_df = do.call("rbind", dawn_chorus_list)

# save csv
write.csv(all_site_avg_df, file = paste0(results_dir, "/averages/site_avg_acoustic_indices.csv"), row.names = FALSE)
write.csv(all_site_by_hour_df, file = paste0(results_dir, "/averages/site_by_hour_acoustic_indices.csv"), row.names = FALSE)
write.csv(all_site_by_hour_all_df, file = paste0(results_dir, "/averages/site_by_hour_all_acoustic_indices.csv"), row.names = FALSE)
write.csv(all_site_dawn_chorus_df, file = paste0(results_dir, "/averages/site_acoustic_indices_dawn_4am-12pm.csv"), row.names = FALSE)
