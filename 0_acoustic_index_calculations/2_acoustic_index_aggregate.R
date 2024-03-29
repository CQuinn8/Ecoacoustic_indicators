# Gather acoustic index results and summarize by site on HPC
# created 27Oct2020
# Colin Quinn

library(data.table)
library(tidyverse)

out_dir = '/results/acoustic_indices_aggregation/by_site/' # directory that will hold by-site csvs

# directory with acoustic index csv results
ai_dir = "/results/acoustic_indices_wavs/" # location of by-wav acoustic index csvs from calculation script
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