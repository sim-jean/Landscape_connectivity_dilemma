"""
Program for Mouysset, Jean, 2023 - Landscape Connectivity Dilemma

STEP 3 : GENERATE MATCHING KEYS FOR EACH BUDGET

Requirement :
- budget in params.py
- optimal successions in results/budget_(params.budget)/

Description :
- When optimal succession is generated, it is not generally part of the unique landscapes we have generated.
- This step generates the matching keys to the unique landscapes generated in STEP 1. It needs to be repeated for each budget specified in the
params file.

Outputs :
- datasets of matching keys for each landscape in results/keys_successions/
"""

import time
import multiprocessing as mp

import modules.params as params
import modules.utilities as utilities
import modules.matcher as matcher

__spec__ = None
__name__ = "__main__"

n_process = 14 #number of processors for multiprocessing

if __name__ == "__main__":
    start = time.time()

    # Apply matching function to whole files
    with mp.Pool(processes=n_process) as pool:
        results = pool.map(matcher.match_heavy_3, list(range(len(params.files))))
    pool.terminate()

    print("The process took " + str((time.time()-start)//60) + ' minutes and ' + str((time.time()-start) % 60) + ' seconds')

