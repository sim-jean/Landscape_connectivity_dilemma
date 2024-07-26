#!/bin/bash
#SBATCH --job-name='stats_4' ### -J 'testJob'
#SBATCH --ntasks=14	            ### -n 1
#SBATCH -p batch	            ### Partition to submit job to 
#SBATCH -o outLog4
#SBATCH -e errLog4
#SBATCH -t 150:00:00

#SBATCH --mail-user=simonjean@ucsb.edu
#SBATCH --mail-type ALL

export PATH=~/anaconda3/bin:$PATH
#####export PATH=/home/zichanghe/anaconda3/envs/qaoa/bin:$PATH

cd $SLURM_SUBMIT_DIR

python statistics_test4.py > outLog4
