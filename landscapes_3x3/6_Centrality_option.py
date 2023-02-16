"""
Program for Mouysset, Jean, 2023 - Landscape Connectivity Dilemma

STEP 6 : DIFFERENT SUCCESSIONS USING EIGENVECTOR CENTRALITY OF GRAPH & STATISTICS

Requirements :
- set budget to value in params.py
- Initial unique landscapes

Description :
Implements succession using a treatment algorithm based on eigenvector centrality of graph. Then computes statistics.

Outputs : datasets of optimal successions & statistics
"""
import pandas as pd

import modules.params as params
import modules.outputs as outputs
import modules.dynamic_programming as dynprog
import multiprocessing as mp
import os
import numpy as np
import networkx as nx


#Load all potential data
files = os.listdir(params.path_to_save_data)

# Loop over files
for k in range(len(params.path_to_save_data)):
    # Load data from file
    lands = dynprog.data_list3(k)
    # Change to np.array
    lands = [np.array(x) for x in lands]

    #Generate succession over 20 periods and array of biodiversity values
    check = dynprog.land_succ(lands)

    # Tweak name for saving
    name = files[k]
    name = name.replace('.pkl','.csv')
    check.to_csv(params.path_to_save_results + "/budget_" + str(params.budget) + "/successions_centrality/succ_cent_budget_"
                 + str(params.budget) + '_' +name)


# Load files resulting from centrality algorithm
files_centrality = os.listdir(params.path_to_save_results+ '/budget_' + str(params.budget) + '/successions_centrality/')

# Compute statistics
for data_path in files_centrality:
    outputs.statistics_dataframe_centrality(data_path)