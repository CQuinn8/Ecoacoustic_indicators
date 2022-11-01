# Purpose : (1) Get unique S2L site names and (2) aggregate CNN inference
# Relies on data calculated and logic in 10.1016/j.ecolind.2022.108831
ll = '/projects/tropics/users/cquinn/R_400/'
.libPaths(c( .libPaths(), ll))
library(data.table)
library(dplyr)


f = "075" # threshold level 
results_dir = '/projects/tropics/users/cquinn/s2l/paper-AcousticIndices/results/ABGQI_inference/'
inference_dir = '/projects/tropics/users/cquinn/s2l/paper-ABGQI_classification/results/cnn_inference_nonXC/'
inference_csvs = list.files(paste0(inference_dir), pattern = "*.csv", full.names = FALSE)

# retrieve all unique site names
rec_name = lapply(inference_csvs, function(x) strsplit(x, "_")[[1]][1]) # get recorder tag
date = lapply(inference_csvs, function(x) strsplit(x, "_")[[1]][2])     # get YYMMDD initial recording
date = lapply(date, function(x) strsplit(x, "-")[[1]][1])               # one site (s2lam057_190619-) is mis-labeled with "-" instead of "_" 
site_names = unique(paste0(rec_name,"_",date))

# get list of files used in ABG CNN training/Val/Test
abg_df = read.csv('/projects/tropics/users/cquinn/s2l/paper-ABGQI_classification/data/cnn_training/verfied_ROIs_thru_5460.csv')
abg_wavs = unique(abg_df$File_name)

# create thresholded class values for ABGQI*
labels = c("Anthro","Bio","Geo","Quiet","Interference")
label_threshes = list()
for(target_label in labels){
  print(target_label)
  # read in optimal threshold based on model optimization
  th_df = read.csv(paste0('/scratch/cq73/projects/S2L/abg_cnn/results/CNN_inference/local_PC_results/IMGNET_S2L_withFreesound-local/performance_fscore_', 
                          f, "/", 
                          target_label, 
                          "_sigmoid_model_acc_metrics.csv"))
  th = subset(x = th_df, X == "th")$test
  
  label_threshes[target_label] = th
}


# iterate through every site's wavs (using only non-overlapping melspec preds (10March21))
error_log <- paste0('/projects/tropics/users/cquinn/s2l/paper-AcousticIndices/results/ABGQI_inference/error_logs.txt')
for(s in 1:length(site_names)){
  temp_site = site_names[s]
  temp_out_path =  paste0(results_dir,"by_site/",temp_site,".csv")
  
  # check if binariaztion has been done
  if(file.exists(temp_out_path)){
    print(paste0(s," : ",temp_site, " has already been binarized... skipping"))
    
  } else  {
    print(paste0(s, ":", temp_site))
    # site specific wavs
    site_wavs = inference_csvs[inference_csvs %like% temp_site]
    
    # check and remove wavs in ABG training step (n = 2367 total)
    site_wavs = site_wavs[!site_wavs %in% abg_wavs]
    
    # iterate through each wav to threshold probabilities
    site_preds = list()
    for(i in 1:length(site_wavs)){
      tryCatch({ 
        temp_csv = read.csv(paste0(inference_dir,site_wavs[i]))
        
        # select only non-overlapping preds for now (10March21)!!
        r = 1:nrow(temp_csv)
        temp_csv = temp_csv[r%%2==1,]
        temp_csv$wav = site_wavs[i] # store wav name 
        
        # binarize each label predictions
        for(j in 1:length(names(label_threshes))){
          temp_th = label_threshes[[j]]
          temp_label = paste0(names(label_threshes)[j],"_bin")
          temp_csv[temp_label] = ifelse(temp_csv[,(j+1)] >= temp_th, 1, 0) # calculate classes based on optimized threshold
        }
        
        temp_csv$site = temp_site
        temp_csv = temp_csv %>%
          select(site, wav, mfcc, Anthro_bin, Bio_bin, Geo_bin, Quiet_bin, Other_bin) %>%
          rename(Anthropophony = Anthro_bin,
                 Biophony = Bio_bin,
                 Geophony = Geo_bin,
                 Quiet = Quiet_bin,
                 Interference = Other_bin)
        
        # store new df in list  
        site_preds[[i]] = temp_csv
        
      }, error = function(e){
        cat("ERROR :",conditionMessage(e), "\n")
        write(toString(temp_site), error_log, append=TRUE)
      })
    }
    # convert list to df
    site_df = do.call("rbind", site_preds)
    write.csv(site_df, row.names = FALSE, temp_out_path)
  }
}