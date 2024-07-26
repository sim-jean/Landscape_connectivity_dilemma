#!/bin/bash
#SBATCH --job-name='5peat4' ### -J 'testJob'
#SBATCH --ntasks=30	            ### -n 1
#SBATCH -p batch	            ### Partition to submit job to 
#SBATCH -o outLog_four
#SBATCH -e errLog_four
#SBATCH -t 150:00:00

#SBATCH --mail-user=simonjean@ucsb.edu
#SBATCH --mail-type ALL

source /home/simonjean/anaconda3/bin/activate r-connectivity

cd $SLURM_SUBMIT_DIR

Rscript Dyn_prog_genetic_server4.R > outLog_four

conda deactivate
