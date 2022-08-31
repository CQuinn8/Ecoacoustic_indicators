# Acoustic Index calculations
# Created: 1/29/19
# Updated: 30/Sept/2019
# Author: Colin Quinn
# Project: Soundscapes 2 Landscapes
# NAU - SICCS

##########################################
# To Do:
## -Dawn chorus timing/elevation needs to be fine-tuned (2/27/2019)
## -create spectrogram output for each site
##    -figure out how to plot HH:MM on ggplot or other method (raster, plot, blah)
## -maybe adjust ggsave sizes now with 7 graphs instead of 6
##########################################

##########################################
## Use this section on monsoon
# # Load libraries
ll = '/home/cq73/R/4_audio/'
library(seewave, lib.loc=ll)
library(soundecology, lib.loc=ll)
library(tuneR, lib.loc=ll)
library(crayon, lib.loc=ll)
library(dplyr, lib.loc=ll)
library(tidyr, lib.loc=ll)
library(reshape2, lib.loc=ll)

# Bring in SLURM parameters
args <-commandArgs(TRUE)
slurm_id <- args[1]
slurm_id <- as.integer(slurm_id) #index for job array
chunk_size <- args[2] # how many files
chunk_size <- as.integer(chunk_size)
wd_out = args[3] # results location
filenames_csv = args[4] # csv listing toDo wavs

print(paste("Slurm ID:", slurm_id))
print(paste("n_files =", chunk_size))
print(paste("Results dir:", wd_out))
print(paste("CSV with ToDo files = ",filenames_csv))

# start and stop files based on slurm ID and chunk size from bash script
start_index = (slurm_id-1) * chunk_size + 1
end_index = slurm_id * chunk_size 
print(paste("Start pos = ", start_index))
print(paste("End pos =", end_index))

# use csv to get list of filenames to do
# need to search for .wav and .WAV b/c both exist
#file_name <- list.files(path = wd, pattern = ".wav|.WAV", full.names = FALSE)
#filenames_csv = "//nau.froot.nau.edu/cirrus/scratch/cq73/projects/S2L/acoustic_indices/results/toDo_2020.csv"
file_names = read.csv(filenames_csv)$wavs
file_list = as.vector(file_names[start_index:end_index]) # slice only wavs in this slurm array task

# list files in site folder
n = length(file_list)

print("\nLoaded audio files list...")
print(paste0("Length of audio file list:", n))
print(paste0("First file = ",start_index,":",file_list[1]))
print(paste0("Last file = ",end_index,":",file_list[n]))


##########################################
# ACOUSTIC INDEX CALCULATIONS #
# Read in wav files and string soundecology and seewave metrics with parameters influenced by Eldridge et al. 2018
  # simply add a new list object below, to index calculation, and to final list index for new acoustic index
  # - currently calculating (ACI,NDSI,BI,AEI,ADI,H,Ht)
  #  -this combination of AI was chosen based on Eldridge et al. 2018, Mammides et al. 2017, Moreno-Gomez et al. 2019, Sueur 2018
  
index_list <- list()
error_log <- paste0("/projects/tropics/users/cquinn/s2l/code/paper1-SpatialMapping/acoustic_indices/error_log_toDo.txt")
# iterate through each sound file
for(i in 1:n){
	temp_wav = file_list[i]
	temp_wav_sans_ext = tools::file_path_sans_ext(temp_wav) #remove extension	
	temp_out_wav = strsplit(temp_wav_sans_ext, '/')[[1]] # split file path
	temp_out_wav = tail(temp_out_wav, n = 1) # file name only
	out_file = paste0(wd_out, temp_out_wav,'.csv') # csv output
	
	# check if wav has been computed
	if (file.exists(out_file)){  
		#output statement and skip 
		cat(paste('File exists: skipping calculations'))
	
	}else{  	  
		tryCatch({ 
		
		  #feed back current file, in case of partial read through
		  cat(paste0(i,"- Working on:", temp_out_wav))
		  
		  # read in wav files from wd
		  soundfile = readWave(temp_wav)
		  samp_freq = soundfile@samp.rate # 44100 for LG; 48000 for AM
		  
		  # calculate each acoustic index 
		  # Soundecology derived indices
		  ndsi_indices <- soundecology::ndsi(soundfile, fft_w = 1024, anthro_min = 1000, anthro_max = 2000, bio_min = 2000, bio_max = 10000)
		  
		  ndsi <- ndsi_indices$ndsi_left
		  
		  anthro <- ndsi_indices$anthrophony_left
		  
		  bio <- ndsi_indices$biophony_left
		  
		  # Biophony Index (Boelman et al., 2007: Appendix A <https://esapubs.org/archive/appl/A017/086/appendix-A.htm>)
		  bi <- soundecology::bioacoustic_index(soundfile, min_freq = 2000, max_freq = 10000)$left_area 
		  
		  aei <- soundecology::acoustic_evenness(soundfile, max_freq = 10000, db_threshold = -50, freq_step = 1000)$aei_left
		  
		  adi <- soundecology::acoustic_diversity(soundfile, max_freq = 10000, db_threshold = -50, freq_step = 1000, shannon = TRUE)$adi_left

		  
	    # SEEWAVE derived indices
		  spec = spec(soundfile, f = samp_freq, plot = FALSE)
		  
		  env <- seewave::env(soundfile, plot = FALSE)
		  
		  aci <- seewave::ACI(soundfile, flim = c(1,10))
		  
		  h <- seewave::H(soundfile, wl = 512)
		  
		  # Temporal entropy (Sueur et al., 2008)
		  ht <- seewave::th(env)
		  
		  # Shannon's spectral entorpy (Sueur et al., 2008)
		  sh <- seewave::sh(spec)
		  
		  rms <- seewave::rms(env)
		  
		  zcr <- seewave::zcr(soundfile, plot = FALSE)
		  
		  M <- seewave::M(soundfile)
		  
		  rough <- seewave::roughness(meanspec(soundfile, plot = FALSE)[,2])
		  
		  sfm <- seewave::sfm(spec)
		  
		  rugo <- seewave::rugo(soundfile@left/max(soundfile@left))
		  
		  ##########################################
		  #Calculate datestamp inputs for dataframe
		  sitename = strsplit(temp_out_wav,'_')[[1]][1]
      YYMMDD = strsplit(temp_out_wav,'_')[[1]][2]
		  YY = substr(YYMMDD,1,2)
		  YYYY = paste0('20',YY)		  
		  MM = substr(YYMMDD,3,4)
		  DD = substr(YYMMDD,5,6)
		  hh_mm = strsplit(temp_out_wav,'_')[[1]][4]
		  hh = substr(hh_mm,1,2)
		  mm = substr(hh_mm,4,5)
		  DDhh = paste0(DD,hh)
		  hhmm = paste0(hh,mm)
		  
		  ##########################################
		  # organize derived products for each file in loop
		  temp_df = data.frame("path" = file_list[i], 
							  'file' = temp_out_wav,
							  "ACI" = aci,
							  "ADI" = adi,
							  "AEI" = aei,
							  'NDSI' = ndsi,
							  'NDSI_A' = anthro,
							  'NDSI_B' = bio,
							  'BI' = bi,
							  'H' = h,
							  'Ht' = ht,
							  'Hs' = sh, # signal noisiness/pureness (similar to sfm)
							  'M' = M,
							  'R' = rough,
							  'sfm' = sfm, # spectral flatness measure (noisy -> 1; pure tone -> 0)
							  'rugo' = rugo, # noisy is higher
							  'zcr_max' = max(zcr[,2]),
							  'zcr_mean' = mean(zcr[,2]),
							  'zcr_min' = min(zcr[,2]),
							  'YYYY' = YYYY, 
							  'MM' = MM, 
							  'DD' = DD, 
							  'hh' = hh, 
							  'mm' = mm, 
							  'DDhh' = DDhh, 
							  'hhmm' = hhmm)
      print(temp_df)     
                   
			# Data Output #########################################
	    # save the slurm task csv
	    print(paste("SAVING OUTPUT FILE:", out_file))
	    write.csv(temp_df, file = out_file, row.names = FALSE)  
         
		}, error = function(e){
		  cat("ERROR :",conditionMessage(e), "\n")
		  write(toString(file_list[i]), error_log, append=TRUE)
		  })  
	}
}

