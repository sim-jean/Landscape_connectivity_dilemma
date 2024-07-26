#!/bin/bash
#SBATCH --job-name='sample' ### -J 'testJob'
#SBATCH --ntasks=0	            ### -n 1
#SBATCH -p batch	            ### Partition to submit job to 
#SBATCH -o outLog_sample
#SBATCH -e errLog_sample
#SBATCH -t 150:00:00

#SBATCH --mail-user=simonjean@ucsb.edu
#SBATCH --mail-type ALL

source /home/simonjean/anaconda3/bin/activate r-connectivity

cd $SLURM_SUBMIT_DIR

Rscript sample_analysis_server.R > outLog_sample

conda deactivate
