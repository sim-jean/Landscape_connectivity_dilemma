#!/bin/bash
#SBATCH --job-name='5peat3' ### -J 'testJob'
#SBATCH --ntasks=30	            ### -n 1
#SBATCH -p batch	            ### Partition to submit job to 
#SBATCH -o outLog_three
#SBATCH -e errLog_three
#SBATCH -t 150:00:00

#SBATCH --mail-user=simonjean@ucsb.edu
#SBATCH --mail-type ALL

source /home/simonjean/anaconda3/bin/activate r-connectivity

cd $SLURM_SUBMIT_DIR

Rscript Dyn_prog_genetic_server3.R > outLog_three

conda deactivate
