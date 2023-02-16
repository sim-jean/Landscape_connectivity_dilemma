"""
Program for Mouysset, Jean, 2023 - Landscape Connectivity Dilemma

STEP 4 : MATCHING SUCCESSIONS WITH VARIOUS BIODIVERSITY LEVELS FOR GIVEN BUDGET

Requirements :
- set budget to value in params.py
- Matching keys in 'keys_succession'

Description:
- For each landscape, successions for 20 periods, with given biodiversity constraints

Outputs : All successions from initial landscapes for 20 periods with all biodiversity constraint levels.

"""
import time
import multiprocessing as mp

import modules.successions as successions
import modules.params as params

# Parameters :
n_processes = 14 #number of processors for multiprocessing

step = 2 #step for biodiversity constraint increment

__name__ = '__main__'
__spec__ = None
if __name__ == "__main__":

    start = time.time()
    #Loop over all files
    for j in params.files:
        # Form candidates with file name and biodiversity level
        candidates = list(zip([j]*len(list(range(2, max(params.biodiv3), step))), list(range(2, max(params.biodiv3), step))))
        # Multiprocessing using starmap (tuple of inputs) and candidates
        with mp.Pool(processes=n_processes) as pool:
            results = pool.starmap(successions.succession_dataframe3, candidates)
        pool.terminate()
        #Following up:
        print('File '+str(j)+' took '+str(time.time()-start)+' seconds')
    print("The process took " + str((time.time()-start)//60) + ' minutes and ' + str((time.time()-start) % 60) + ' seconds')