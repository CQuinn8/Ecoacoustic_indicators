
################## PAIRED ARU FUNCTIONS #######################
# read in paired ABGQI data
# read in all pairs and create dataframes of each paired set of predictions by overlapping time
abgqi_pair_reader = function(paired_df){
  pair_list = list()
  for(pair in seq_along(paired_df$SiteID)){
    temp_row = paired_df[pair,]
    print(paste0("Pair: ",pair))
    
    #read in csvs
    am = read.csv(paste0(wd, 'results/paired_ARUs/paired_ABGQI/', temp_row$SiteID, '.csv')) %>%
      dplyr::select(Anthro, Bio, Geo, Other, Quiet, wav, mfcc) %>%
      dplyr::rename(Anthrophony = Anthro, Biophony = Bio, Geophony = Geo, Interference = Other, Quiet = Quiet, wavAM = wav) %>%
      gather(variable, value, -wavAM, -mfcc)
    lg = read.csv(paste0(wd, 'results/paired_ARUs/paired_ABGQI/', temp_row$Paired_site, '.csv')) %>%
      dplyr::select(Anthro, Bio, Geo, Other, Quiet, wav, mfcc) %>%
      dplyr::rename(Anthrophony = Anthro, Biophony = Bio, Geophony = Geo, Interference = Other, Quiet = Quiet, wavLG = wav) %>%
      gather(variable, value, -wavLG, -mfcc)
    
    # build dateframe
    am$site_AM = temp_row$SiteID
    am$DD = substr(am$wav, 25, 26)
    am$HH = substr(am$wav, 28, 29)
    am$mm = substr(am$wav, 31, 32)
    am$DDHHmm = paste0(am$DD, am$HH, am$mm)
    
    lg$site_LG = temp_row$Paired_site
    lg$DD = substr(lg$wav, 25, 26)
    lg$HH = substr(lg$wav, 28, 29)
    lg$mm = substr(lg$wav, 31, 32)
    lg$DDHHmm = paste0(lg$DD, lg$HH, lg$mm)
    
    pair_list[[pair]] = am %>%
      inner_join(lg, by = c('DD','HH','mm','variable','mfcc')) %>%
      dplyr::select(wavAM, mfcc, variable, value.x, value.y, DDHHmm.x, site_AM, site_LG) %>%
      dplyr::rename(wav = wavAM, melspec = mfcc, soundType = variable, DDHHmm = DDHHmm.x, AM = value.x, LG = value.y) %>%
      mutate(SiteID = pair)
    
  }
  pair_preds = do.call(rbind, pair_list)
}

# read in acoustic index data
  acoustic_index_pair_reader = function(paired_df){
  pair_list = list()
  for(pair in seq_along(paired_df$SiteID)){
    temp_row = paired_df[pair,]
    print(paste0("Pair: ",pair))
    
    #read in csvs
    am = read.csv(paste0(wd, 'results/paired_ARUs/paired_acoustic_indices/', temp_row$SiteID, '.csv')) %>%
      select(-YYYY, -MM, -DD, -hh, -mm, -DDhh, -hhmm) %>%
      rename(site_AM = site, wav = file) %>%
      gather(variable, value, -wav, -site_AM)
    lg = read.csv(paste0(wd, 'results/paired_ARUs/paired_acoustic_indices/', temp_row$Paired_site, '.csv')) %>%
      select(-YYYY, -MM, -DD, -hh, -mm, -DDhh, -hhmm) %>%
      rename(site_LG = site, wav = file) %>%
      gather(variable, value, -wav, -site_LG)
    
    # build dateframe
    am$DD = substr(am$wav, 25, 26)
    am$HH = substr(am$wav, 28, 29)
    am$mm = substr(am$wav, 31, 32)
    am$DDHHmm = paste0(am$DD, am$HH, am$mm)
    
    lg$DD = substr(lg$wav, 25, 26)
    lg$HH = substr(lg$wav, 28, 29)
    lg$mm = substr(lg$wav, 31, 32)
    lg$DDHHmm = paste0(lg$DD, lg$HH, lg$mm)
    
    pair_list[[pair]] = am %>%
      inner_join(lg, by = c('DD','HH','mm','variable')) %>%
      dplyr::select(wav.x, variable, value.x, value.y, DDHHmm.x, site_AM, site_LG) %>%
      dplyr::rename(wav = wav.x, acIndex = variable, DDHHmm = DDHHmm.x, AM = value.x, LG = value.y) %>%
      mutate(SiteID = pair)
    
  }
  pair_indices = do.call(rbind, pair_list)
}


# boxplots for ARU comparison
index_boxplot = function(df){
  df %>%
    ungroup() %>%
    select(AM, LG) %>%
    gather(ARU, value) %>%
    ggplot(aes(x = ARU, y = value)) + 
    geom_jitter(width = 0.1, alpha = 0.1) + 
    geom_boxplot(alpha = 0.4, outlier.shape = NA, notch = TRUE)
}

# corr plots for ARU comparison
corr_plotter = function(df){
  df%>%
    ungroup() %>%
    select(where(is.numeric)) %>%
    select(contains(c("AM","LG"))) %>%
    GGally::ggpairs(aes(alpha = 0.05), progress = FALSE)
}

# function to normalize data from 0 to 1
min_max_norm <- function(x, min_x, max_x) {
  (x - min_x) / (max_x - min_x)
}

center_scale = function(x , mu, sd){
  (x - mu) / sd
}


# Mean independent (ie robust to outlier) estimate of scale for an array of values
MAD_z = function(array){
  med = median(array)
  n = length(array)
  
  # derive median of all absolute deviations from the median
  MAD_med = median(sapply(array, function(x) abs(x - med)))
  
  # apply MAD correction for unbiased normal distribution
  MAD = 1.483 * MAD_med
  
  # median and MAD estimator of robust scores
  est = sapply(array, function(x) (x - med)/MAD)
}

# IDentify outlier values using the MAD approach and threshold score at 99 percentile
outlier_id = function(array){
  MAD_values = MAD_z(array)
  threshold = quantile(abs(MAD_values), probs = 0.99)
  outlier_positions = which(abs(MAD_values) >= threshold)
  return(outlier_positions)
}

toLogit <- function(x){
  cx<-ifelse(x==0,0.00001,
             ifelse(x==1,0.99999,x))
  lgx<-log(cx /(1-cx))
  return(lgx)
}

# x is the value to logit-transform from
fromLogit<-function(x){
  bt<-exp(x)/(1+exp(x))
  return(bt)
}


# GAM prediction visualization
pred_plot = function(data_df, model_fit, index_name, link) {
  # predict based on original values
  data_df$y_pred = predict(model_fit, newdata = data_df, type = "response")
  
  # covariates in model
  vars = str_replace(names(model_fit$sp), "s*", "")
  vars = gsub("[()]", "", vars)
  
  data_df %>%
    select(vars, contains(index_name), y_pred) %>%
    rename(acoustic_index = contains(index_name, ignore.case = FALSE)) %>%
    gather(covariate, value, -acoustic_index, -y_pred, -ARU) %>%
    ggplot(aes(x = value, y = y_pred)) +
      geom_point(alpha = 0.4) +
      geom_point(aes(x = value, y = acoustic_index, colour = ARU), alpha = 0.3) +
      facet_wrap(~covariate) +
      ggtitle(index_name)
  
  # data_df %>%
  #   select(-site,-wavs) %>%
  #   rename(acoustic_index = contains(index_name)) %>%
  #   mutate(y_diff = acoustic_index - y_pred) %>%
  #   gather(covariate, value, -acoustic_index, -y_pred, -y_diff, -ARU) %>%
  #     ggplot(aes(x = value, y = y_diff, alpha = 0.4)) +
  #     geom_point() +
  #     facet_wrap(~covariate, scales = 'free_x')
}

# GAMs
# first derivative (i.e., partial effects slopes)
slope_summary = function(model) {
  # covariates in model
  vars = data.frame(vars = names(model$sp)) %>%
    filter(!grepl("ARU", vars, ignore.case = TRUE))
  
  # derivatives
  fd = lapply(vars$vars, function(x) derivatives(model, term = x, partial_match = TRUE))
  fd = do.call("rbind", fd)
  
  return(fd)
}
