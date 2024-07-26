#!/bin/bash

#SBATCH --job-name='Py_MultiPro'    ### -J 'testJob'

#SBATCH --ntasks=10	            ### how many cores you request

#SBATCH -p batch	            ### Partition to submit job to 

#SBATCH -e errLog                   ### generate an error file (not necessary, but recommend)
#SBATCH -t 00:10:00



#SBATCH --mail-user=simonjean@ucsb.edu   ### Send you a reminder email (not necessary, but recommend)
#SBATCH --mail-type ALL                   ### when the job start and end (not necessary, but recommend)



##### export your python environment

export PATH=/home/simonjean/software/anaconda3/bin:$PATH





cd $SLURM_SUBMIT_DIR



python test_simon.py > outLog

