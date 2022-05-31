ll = '/projects/tropics/users/cquinn/R_400/'
.libPaths(c( .libPaths(), ll))
library(dplyr)

wd = '/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/'

# read in paired ARU csv with AM site list and LG site list
paired_df = read.csv(paste0(wd,'results/paired_ARUs/colocated_ARUs-paired.csv')) %>%
  dplyr::select(SiteID, Number.of.Recordings, Paired_site) %>%
  filter(grepl('am', SiteID)) # 24 paired sites (48 ARUs)

acoustic_index_wd = paste0(wd, 'results/acoustic_indices_aggregation/by_site/')

# Get list of acoustic index csvs by site
paired_index_csvs = paste0(acoustic_index_wd, c(paired_df$SiteID, paired_df$Paired_site), '.csv')

# copy csvs to paired aru dir
to_dir = paste0(wd, 'results/paired_ARUs/paired_acoustic_indices/')
file.copy(from = paired_index_csvs, to = to_dir)
