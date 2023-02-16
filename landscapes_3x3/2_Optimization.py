"""
Program for Mouysset, Jean, 2023 - Landscape Connectivity Dilemma

STEP 2 : OPTIMIZATION GIVEN BUDGET

Requirements :
- set budget to value in params.py

Description :
Script to optimize the succession of unique landscapes:
- For each landscape, optimal successions are derived
- For each biodiversity level, optimal successions are derived.

Outputs : datasets of optimal successions with unmatched output
"""

import time
import multiprocessing as mp

import modules.utilities as utilities
import modules.dynamic_programming as dynprog
import modules.params as params


import os
import re
import pandas as pd
__name__ = "__main__"
__spec__ = None

n_process = 14 #number of processors for multiprocessing

if __name__ == "__main__":
    # Get files on which to perform optimization from path to save data
    files = os.listdir(params.path_to_save_data)

    debut2 = time.time()
    print("Process is running:")
    # Set names for data storage
    list_names = [["land"], ["success_biod_" + str(i) for i in params.biodiv3], ["fire_biod_" + str(i) for i in params.biodiv3],
                  ["value_biod_" + str(i) for i in params.biodiv3]]
    flat_list_names = [item for sublist in list_names for item in sublist]

    # Loop over files:
    for i in list(range(len(files))):
        size = 4
        #load data environment using data_list
        list_nz = dynprog.data_list3(i)

        start = time.time()
        with mp.Pool(processes=n_process) as pool:
            results = pool.map(dynprog.lowlev3_dynprog,
                                   list_nz)
        pool.terminate()

        #region save results
        d = pd.DataFrame(results, columns=list(flat_list_names))

        # Opt succession
        indexes = re.findall(r'\d+', str(files[i]))

        d.to_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(params.budget) + "/land3_budget_" + str(params.budget) + '_' + str(
                indexes[1]) + ".csv", index=False)
        print(str(time.time() - start) + " seconds elapsed for " + str(indexes[1]) + " non-zeros ")
        #endregion
    print("Overall :" + str(time.time() - debut2) + "seconds have elapsed")

