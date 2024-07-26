#!/bin/bash
#SBATCH --job-name='stats_1' ### -J 'testJob'
#SBATCH --ntasks=15	            ### -n 1
#SBATCH -p batch	            ### Partition to submit job to 
#SBATCH -o outLogf
#SBATCH -e errLogf
#SBATCH -t 150:00:00

#SBATCH --mail-user=simonjean@ucsb.edu
#SBATCH --mail-type ALL

export PATH=~/anaconda3/bin:$PATH
#####export PATH=/home/zichanghe/anaconda3/envs/qaoa/bin:$PATH

cd $SLURM_SUBMIT_DIR

python statistics_test5.py > outLogf
