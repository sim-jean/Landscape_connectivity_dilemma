#!/bin/bash
#SBATCH --job-name='lf' 
#SBATCH --ntasks= 40	            ### -n 1
#SBATCH -p batch	            ### Partition to submit job to 
#SBATCH -o outLoga
#SBATCH -e errLoga
#SBATCH -t 150:00:00

#SBATCH --mail-user=simonjean@ucsb.edu
#SBATCH --mail-type ALL

source /home/simonjean/anaconda3/bin/activate r-connectivity

cd $SLURM_SUBMIT_DIR

Rscript large_landscapes_analysis_5p.R > outLog5

conda deactivate
