import numpy as np
import os
import time
import pandas as pd
import params_server4 as params
import outputs_server4 as outputs
import utilities_server as utilities
import math
import multiprocessing as mp
import csv

n_cores = os.getenv('SLURM_NTASKS','1')
n_cores = int(n_cores)
print(str(n_cores) + " cores are running with budget :" + str(params.budget))
print("Starting time :", str(time.ctime(time.time()) ))

step = 10
__name__ = '__main__'
__spec__ = None

if __name__ == "__main__":
    start = time.time()
    candidates = os.listdir("/home/simonjean/data/budget_" + str(params.budget) + '/statistics')

#    file = open('/home/simonjean/data/budget_' + str(params.budget) +'/update.csv', 'w+', newline='')
    #writing the data into the file
#    with file as f:
#        write = csv.writer(f)
#        write.writerows(candidates)

#    with open('/home/simonjean/data/budget_' + str(params.budget) + '/update.csv', 'r') as my_file:
#        reader = csv.reader(my_file, delimiter = '\t')
#        candidates = list(reader)
#    candidates = [x[0].replace(',','') for x in candidates]
    #candidates = candidates[0:1512]
    #candidates = candidates[1512:3024]
    #candidates = candidates[3024:4536]
    #candidates = candidates[4536:]
    print(len(candidates))
    #print(candidates)
    with mp.Pool(processes = n_cores) as pool:
        results = pool.map(outputs.cv_time_period_by_cycle, candidates)
    print("Process is over and took :" +str(time.time()-start) + "seconds" )
    #for candidate in candidates:
    #   outputs.statistics_dataframe(candidate)
    #   print(candidate + " is done!")
