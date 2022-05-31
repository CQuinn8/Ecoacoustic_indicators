# Gather acoustic index resultsand summarize by site on HPC
# created 27Oct2020
# Colin Quinn

ll = '/projects/tropics/users/cquinn/R_400/'
library(crayon, lib.loc = ll)
library(dplyr, lib.loc = ll)
library(data.table)
library(tidyverse, lib.loc = ll)

out_dir = '/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/acoustic_indices_aggregation/by_site/'

# directory with acin results
ai_dir = "/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/acoustic_indices_wavs/"
ai_files = list.files(ai_dir, pattern = "\\.csv$")
ai_files = data.frame(wavs = ai_files)

# get unique site names from files
temp_rec = sapply(ai_files$wavs, function(x) strsplit(x, '_')[[1]][1])
temp_sites = sapply(ai_files$wavs, function(x) strsplit(x, '_')[[1]][2])
temp_sites = paste0(temp_rec, '_', temp_sites)

# clean up any sites with errors
unique_sites = sapply(temp_sites, function(x) strsplit(x, '-')[[1]][1])
unique_sites = unique(unique_sites)


# read in all csvs with that site name
for(i in seq_along(unique_sites)){
  temp_site = unique_sites[i]
  temp_out = paste0(out_dir, temp_site, '.csv')
  
  if(file.exists(temp_out) == TRUE){
    print(paste0("ALREADY COMPLETE ", i,': ', temp_site))
  }  else{
    print(paste0("Working on site ", i,': ', temp_site))
    
    # use current site to subset acoustic index wav csvs
    site_wavs = ai_files[grepl(temp_site, ai_files$wavs), , drop = FALSE]$wavs
    site_wavs_path = paste0(ai_dir, site_wavs)
    
    # read in all csvs and concatenate into single dataframe
    tbl_fread = do.call(rbind, lapply(site_wavs_path, fread))
    tbl_fread$site = temp_site
    tbl_fread = tbl_fread[,-1] # remove full path metadata
    
    # save site dataframe
    write.csv(tbl_fread, file = temp_out, row.names = F)
  }
}





# csv with site subset
# site_csv = "/scratch/cq73/projects/S2L/acoustic_indices/data/site_subset.csv"
# site_df = read.csv(site_csv)

# acin_list = list()
# error_files = list()
# index = 1
# for(i in 1:length(ai_files)){
#   if(i %% 10000 == 0){
#     print(paste(i, ":", ai_files[i]))
#   }
#   
#   tryCatch({
#     # read in csv
#     temp = read.csv(paste0(ai_dir, ai_files[i]))
#     acin_list[[i]] = temp
#     
#   }, error = function(e){
#     cat("ERROR :",conditionMessage(e), "\n")
#     error_files[[index]] = ai_files[i]
#     index = index + 1
#   })
# }

# # search list for sites
# site_df = as.character(site_df$V1)
# site_acin = list()
# index = 1
# for(i in 1:length(acin_list)){
#   temp = acin_list[[i]]
#   
#   for(j in 1:length(site_df)){
#     site = site_df[j]    
#     if(any(grepl(site, temp$file, fixed = TRUE)) == TRUE){
#       site_acin[[index]] = temp[grep(site, temp$file),]
#       index = index + 1
#     }
#   }
#   if(i %% 250 == 0){
#     print(i)
#   }
# }
# 
# acin_df = dplyr::bind_rows(site_acin)
# df_sites = vector()
# for(row in 1:nrow(acin_df)){
#   temp = acin_df[row,]
#   temp = unlist(strsplit(as.character(temp$file), "_"))[c(1,2)]
#   df_sites[row] = paste0(temp[1], "_", temp[2])
# }
# 
# acin_df$site = df_sites
# 
# # group by site and get average
# x = acin_df %>%
#   group_by(site)
