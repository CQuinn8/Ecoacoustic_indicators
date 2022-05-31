#!/bin/bash

# acoustic_index.sh
# ----- ABOUT -----
# By Colin Quinn, NAU
# cq73@gmail.com
# Created: 2/28/2019
# Last updated: 30-Nov-2021
#
# Purpose: 
# Calculates acoustic indices (n=15+) for directories of wav files using slurm tasks for n_files
# 
# NOTE: 
#	1) run ls | wc -l in data dir to get total files
#
# TODO:
# 


# ENTER SETUP SCRIPT
# - SLURM -
#SBATCH --time=00:02:00
#SBATCH --output='/scratch/cq73/projects/S2L/acoustic_indices/code/inference.out'

echo "Entering first slurm script"
# --- SLURM Settings ---
arrayLim=2000 # how many array jobs to have at once
files_per_task=1 # 50; how large each job should be, number of files
timeReq="00:03:00" # 18:00
memoryReq="5GB"

# where output, first .sh script will be written to
script_dir='/projects/tropics/users/cquinn/s2l/code/paper1-SpatialMapping/acoustic_indices/'

# --- DATA & RESULT DIR Settings ---
results_dir='/projects/tropics/users/cquinn/s2l/paper1-SpatialMapping/results/acoustic_indices/'
num_files=163
toDo_csv="/projects/tropics/users/cquinn/s2l/code/paper1-SpatialMapping/acoustic_indices/toDo_wavs.csv"

# --- Data preprocessing ---
date_time=`date +%Y%m%d_%H%M%S`
# num_files=$(find $data_dir -type f | wc -l) # the number of files, 
echo "Working data dir : "$data_dir
echo "Number of files = "$num_files

ntask=$(echo "scale=2 ; $num_files / $files_per_task" | bc) # calculate how many jobs there should be
echo "Number of tasks ="$ntask

njobs=$(echo "($ntask+0.999)/1" | bc)
echo "Number of jobs ="$njobs


# Specify the basename for the output log file
slurmLogName=$script_dir'logs/%a.out'


# --- SLURM ARRAY ---
# Create the slurm array job script in the slurm folder
cd $script_dir'sbatch/'
scriptName='sbatch_acoustic_index'$date_time'.sh'

cat > $scriptName <<EOT
#!/bin/bash

# - SLURM -
#SBATCH --job-name=AcIn                  # defined name for jobstats
#SBATCH --output=$slurmLogName           # location for .out file
#SBATCH --partition=core                 # partition name
#SBATCH --time=$timeReq                  # walltime defined above
#SBATCH --cpus-per-task=1
#SBATCH --mem=$memoryReq                 # mem in GB
#SBATCH --array=[1-$njobs]%$arrayLim


start=$(date +%s)
date_time_inner=`date +%Y%m%d_%H%M%S`
echo "The starting date_time: " \$date_time_inner
echo
echo "SLURM_JOBID: "\$SLURM_JOBID
echo "SLURM_ARRAY_JOB_ID: "\$SLURM_ARRAY_JOB_ID
echo "SLURM ARRAY TASK ID: "\$SLURM_ARRAY_TASK_ID
echo

# - MODULES -
module load R

echo
echo "---------Entering Rscript---------"
# - R SCRIPT -
Rscript /projects/tropics/users/cquinn/s2l/code/paper1-SpatialMapping/acoustic_indices/acoustic_index_15_toDo.R \$SLURM_ARRAY_TASK_ID $files_per_task $results_dir $toDo_csv

# - ENDING -
echo "Ended at:"
date
echo
end=`date +%s`
totTime=$((end-start))
echo Total time: $totTime sec

EOT

# Run the slurm array job script
sbatch $scriptName