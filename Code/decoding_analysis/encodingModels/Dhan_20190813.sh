#!/bin/bash
# Dhan_20190813.sbatch 115 sess
#SBATCH -J Dhan13Enc
#SBATCH -p olveczky,serial_requeue,shared,general   # partition (queue)
#SBATCH --exclude=shakgpu[01-50]
#SBATCH -N 1                 # number of nodes
#SBATCH -c 8            # number of cores
#SBATCH --mem 16000          # memory for all cores
#SBATCH -t 6-00:00           # time (D-HH:MM)
#SBATCH --export=ALL
#SBATCH -o Dhan13_%A_%a.%N.%j.out     # STDOUT
#SBATCH -e Dhan13_%A_%a.%N.%j.err     # STDERR
 
module load matlab/R2018a-fasrc01
cd /n/holylfs02/LABS/olveczky_lab/Ashesh/DLS_Analysis/Code/

# Create a local work directory
mkdir -p /scratch/dhawale/$SLURM_JOB_ID

matlab -nodisplay -nosplash -r "fitEncodingModelsRCWrapper('Dhanashri', ${SLURM_ARRAY_TASK_ID}, '../Analysis/Encoding/20190813/')"

# Cleanup local work directory
rm -rf /scratch/dhawale/$SLURM_JOB_ID