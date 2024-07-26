#!/bin/bash
#SBATCH --job-name='3peat Heuristics' ### -J 'testJob'
#SBATCH --ntasks=20	            ### -n 1
#SBATCH -p batch	            ### Partition to submit job to 
#SBATCH -o outLogf
#SBATCH -e errLogf
#SBATCH -t 150:00:00

#SBATCH --mail-user=simonjean@ucsb.edu
#SBATCH --mail-type ALL

source /home/simonjean/anaconda3/bin/activate r-connectivity

cd $SLURM_SUBMIT_DIR

Rscript Dyn_prog_genetic_server.R > outLogf

conda deactivate
