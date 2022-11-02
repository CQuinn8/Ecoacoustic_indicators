# function to normalize data from 0 to 1 based on min and max
min_max_norm <- function(x) {
  min_x = min(x) - 1e-7
  max_x = max(x)
  (x - min_x) / (max_x - min_x)
}

# Calculates mean and sd in tidy summary
mean_sd <- list(
  mean = ~mean(.x, na.rm = TRUE), 
  sd = ~sd(.x, na.rm = TRUE)
)

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

# Function that converts to logit space
toLogit <- function(x){
  cx<-ifelse(x==0,0.00001,
             ifelse(x==1,0.99999,x))
  lgx<-log(cx /(1-cx))
  return(lgx)
}

# Convert to original number line from logit space
# x is the value to logit-transform from
fromLogit<-function(x){
  bt<-exp(x)/(1+exp(x))
  return(bt)
}

# GAM prediction visualization
pred_plot = function(data_df,     # dataframe containing response and covs
                     model_fit,   # GAM fit object
                     index_name,  # char matching response name
                     link) {      # if applicable, any transformation for visual
  
  # predict based on original values
  data_df$y_pred = predict(model_fit, newdata = data_df, type = "response")
  
  # covariates in model
  vars = str_replace(names(model_fit$sp), "s*", "")
  vars = gsub("[()]", "", vars)
  
  data_df %>%
    select(vars, index_name, y_pred, ARU) %>%
    rename(acoustic_index = index_name) %>%
    gather(covariate, value, -acoustic_index, -y_pred, -ARU) %>%
    ggplot(aes(x = value, y = y_pred)) +
      geom_point(alpha = 0.4) +
      geom_point(aes(x = value, y = acoustic_index, colour = ARU), alpha = 0.3) +
      facet_wrap(~covariate, scales = 'free') +
      ggtitle(paste0(index_name, ": black = predicted values"))
}

# GAM first derivative calcualtion (i.e., partial effects slopes)
slope_summary = function(model) {
  # covariates in model
  vars = data.frame(vars = names(model$sp)) %>%
    filter(!grepl("ARU", vars, ignore.case = TRUE))
  
  # derivatives for each covariate 
  fd = lapply(vars$vars, function(x) derivatives(model, term = x, partial_match = TRUE))
  fd = do.call("rbind", fd)
  
  return(fd)
}

# Formated index names used for plotting
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

# helper function that looks up unformatted index name with formatted
index_labeller <- function(variable,value){
  return(index_names[value])
}

# PDP plot function that helps format correct number of rows based on number of covariates
write_pdp <- function(plot_obj, index_name, rows = 1, out_dir){
  
  # custom figure height based on number of covariates divided by 4
  if(rows == 1){
    h = 3
  } else {
    h = 6
  }
  
  # save ggplot object
  ggsave(filename = paste0(index_name, '_pdp.png'), 
         plot = plot_obj, 
         device = 'png',
         path = out_dir,
         width = 8, height = h, dpi = 500)
}

