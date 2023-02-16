"""
Program for Mouysset, Jean, 2023 - Landscape Connectivity Dilemma

STEP 5 : STATISTICS ON LANDSCAPES

Requirements :
- set budget to value in params.py
- Successions matched for various biodiversity constraints

Description :
Script to perform the computation of relevant statistics for analysis. Statistics are described in the statistics.py module

Outputs : datasets of optimal successions and associated statistics
"""
import pandas as pd

import modules.params as params
import modules.outputs as outputs
import multiprocessing as mp
import os


__name__ = "__main__"
__spec__ = None

n_process = 14 # Number of processors for multiprocessing

if __name__ == "__main__":

    #Get files from matched successions
    files = os.listdir(
            "C:/Users/jean/PycharmProjects/connectivity_dilemma/data" + str(params.size) + "/results/budget_" + str(
            params.budget) + "/successions_matched/")

    with mp.Pool(processes=n_process) as pool:
        result = pool.map(outputs.statistics_dataframe,files)
    pool.terminate()
    print('Statistics are done for budget ' + str(params.budget))
