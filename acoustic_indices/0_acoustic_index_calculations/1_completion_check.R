# library(tuneR)
# Build acoustic index to do csvs

wd_pre = '/scratch/cq73/projects/S2L/distributionData/'
wd_post = '/scratch/cq73/projects/S2L/audio/'
results_dir = '/projects/tropics/users/cquinn/s2l/paper1-SpatialMapping/results/acoustic_indices/'

# # check audio files
# temp = readWave(paste0(wd_post, 's2lam020_210516_2021-05-16_12-40.WAV'))
# temp@samp.rate

# build csv with all files
# pre_files = list.files(path = wd_pre, pattern = ".wav|.WAV", full.names = TRUE)
# post_files = list.files(path = wd_post, pattern = ".wav|.WAV", full.names = TRUE)
# all_wavs = c(pre_files, post_files)
# all_wavs = sub('//shares.hpc.nau.edu/cirrus','',all_wavs)
# all_wavs = data.frame('wavs' = all_wavs)
# write.csv(all_wavs, '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/code/paper1-SpatialMapping/acoustic_indices/toDo_wavs.csv', row.names = FALSE)

# check for wav completion
pre_files = list.files(path = wd_pre, pattern = ".wav|.WAV", full.names = FALSE)
post_files = list.files(path = wd_post, pattern = ".wav|.WAV", full.names = FALSE)
wavs = c(pre_files, post_files)
results = list.files(path = results_dir, pattern = '*.csv$', recursive = FALSE, full.names = FALSE)

# remove .wav and .pkl sub('.csv','',results)
results_no_ext = tools::file_path_sans_ext(results)

# check pre 2021 completion
pre_no_ext = tools::file_path_sans_ext(pre_files)
pre_toDo = setdiff(pre_no_ext, results_no_ext)

# check post 2021 completion
post_no_ext = tools::file_path_sans_ext(post_files)
post_toDo = setdiff(post_no_ext, results_no_ext)

# add extensions and pathways
pre_wavs_toDo = paste0(wd_pre, pre_toDo, '.wav')

# need more complex 2020 + 2021 check bc files are .wav and .WAV
post_wavs_toDo = list()
for(i in seq_along(post_toDo)){
  print(i)
  f = post_toDo[i]
  temp_f = post_files[post_files %like% f]
  post_wavs_toDo[[i]] = paste0(wd_post, temp_f)
}
post_wavs_toDo = do.call('rbind', post_wavs_toDo)
length(post_wavs_toDo) == length(post_toDo) # check we got all wavs

# write toDo csv
all_wavs = c(pre_wavs_toDo, post_wavs_toDo)
toDo_wavs = data.frame('wavs' = all_wavs)
write.csv(toDo_wavs, '/projects/tropics/users/cquinn/s2l/code/paper1-SpatialMapping/acoustic_indices/toDo_wavs.csv', row.names = FALSE)
