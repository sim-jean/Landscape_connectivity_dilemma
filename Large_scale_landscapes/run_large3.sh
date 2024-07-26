#!/bin/bash
#SBATCH --job-name='Large3' ### -J 'testJob'
#SBATCH --ntasks=10	            ### -n 1
#SBATCH -p batch	            ### Partition to submit job to 
#SBATCH -o outLog3
#SBATCH -e errLog3
#SBATCH -t 150:00:00

#SBATCH --mail-user=simonjean@ucsb.edu
#SBATCH --mail-type ALL

source /home/simonjean/anaconda3/bin/activate r-connectivity

cd $SLURM_SUBMIT_DIR

Rscript large_landscapes_analysis_3p.R > outLog3

conda deactivate
