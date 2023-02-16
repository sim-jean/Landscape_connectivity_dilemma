'''
Module to implement the matching the results of the optimization procedure to their unique landscape counterparts

Functions :
------------
- nb2, nonzero : utility functions for file location
- match_heavy : matching function for 3x3 landscapes
- match_heavy_3 : matching function for 4x4 landscapes
'''

import modules.utilities as utilities
import modules.params as params
#region Packages
import statistics

import numpy as np
import pandas as pd
import time
import ast
import csv
import statistics as stats
import random as rd
import re as re
import math
import multiprocessing as mp
import os
import json
import pickle
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
#endregion

def nonzero(land):
    return sum(land > 0)

def nb2(land):
    return sum(land == 2)

def match_heavy(k):
    '''
    Matches results to their equivalent landscape in the set of unique landscapes for 4x4 landscapes

    :param k: int
    :return: dataset in .csv
    '''

    # Open file
    file = open(params.path+params.files[k])
    nums = re.findall(r'\d+', params.files[k])
    csvreader = csv.reader(file)
    header = next(csvreader)

    # Load data as list
    rows = []
    for row in csvreader:
        row_add = row[1:42]
        rows.append(row_add)
    landscapes = [item for sublist in rows for item in sublist]

    # keep only unique landscapes
    singular_land = [np.array(ast.literal_eval(i), int) for i in list(set(landscapes))]


    path2 = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data"
    files_data = os.listdir(path2)

    #region Find number of twos and ones : will need to load data from unique landscapes with corresponding number of 1 and 2
    twos = set([nb2(i) for i in singular_land])
    nonzeros = set([nonzero(i) for i in singular_land])

    possible = list(itertools.product(nonzeros, twos))

    store = {}
    for i in range(len(possible)):
        store[possible[i]] = [j for j in singular_land if (nonzero(j) == possible[i][0] and nb2(j) == possible[i][1])]
    store2 = {k: v for k, v in store.items() if v}
    del store
    #endregion

    #region Matching & replacing procedure
    debut = time.time()
    replacement = {}

    # Loop over all potential combinations of 1 and 2 that is needed
    for j in list(store2.keys()):
        # If file exists and is not part of a set of combinations of 1 and 2 that had to be cut for better size management,
        # then load it
        if os.path.exists(path2 + "/land4_" + str(j[0]) + "_" + str(j[1]) + ".pkl"):
            filename = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data/land4_" + str(j[0]) + "_" + str(
                j[1]) + ".pkl", 'rb')
            set_to_check = pickle.load(filename)
            filename.close()

            # For each land with corresponding number of ones and twos (stored in j):
            for v in store2[j]:
                # Generate equivalent landsapes
                set_verif = utilities.equi_landscape(utilities.to_matrix(v))
                set_verif.add(tuple(v))
                # Check if there is an intersection with loaded data
                list_verif = [*set_verif.intersection(set_to_check)]
                # Keep it if exists
                replacement[str(tuple(v))] = str(list_verif[0])

        # If the number of 1 and 2s is linked to a set of lands with 1 and 2s such that it had to be cut, load all of the files
        # and repeat the operation.
        elif os.path.exists(path2 + "/land4_" + str(j[0]) + "_" + str(j[1]) + "cut_0.pkl"):
            set_to_check_merg = []
            to_load = [x for x in files_data if x.startswith('land4_' + str(j[0]) + "_" + str(j[1]))]

            for x in to_load:
                filename = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data/"+x, 'rb')
                set_to_check = pickle.load(filename)
                filename.close()
                set_to_check = list(map(tuple, set_to_check))
                set_to_check_merg.append(set_to_check)

            set_to_check_merg = [item for sublist in set_to_check_merg for item in sublist]
            set_to_check_merg = set(set_to_check_merg)

            for v in store2[j]:
                set_verif = utilities.equi_landscape(utilities.to_matrix(v))
                set_verif.add(tuple(v))
                list_verif = [*set_verif.intersection(set_to_check_merg)]
            # new = list(list_verif[0])

                replacement[str(tuple(v))] = str(list_verif[0])

            del set_to_check
        else:
            print('ERROR : no file')
    #endregion

    # Once every land has been examined, replace it with its replacement value
    for i in list(range(len(rows))):
        rows[i] = [replacement.get(n, n) for n in rows[i]]

    # Put in dataframe
    tester3 = pd.DataFrame(rows, dtype=str)

    # depending on the number of numbers in the data_path, save it accordingly like equiv_results_3_2.csv or equiv_results_8_5_cut5.csv
    if len(nums) > 4:
        tester3.to_csv(params.path + "/keys_succession/equiv_results_" + nums[2] + "_" + nums[3] + "_cut_" + nums[4] + '.csv',
                  index=False)
    else:
        tester3.to_csv(params.path + "/keys_succession/equiv_results_" + nums[2] + "_" + nums[3] + ".csv", index=False)

    return print('Dataset ' + str(params.files[k]) +" took " + str(time.time()-debut) + ' seconds')

def match_heavy_3(k):
    '''
    Matches results to their equivalent landscape in the set of unique landscapes for 3x3 landscapes with given budget

    :param k: int
    :return: dataset in .csv
    '''
    file = open(params.path + params.files[k])
    nums = re.findall(r'\d+', params.files[k])
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows = []
    for row in csvreader:
        row_add = row[0:26]
        rows.append(row_add)
    landscapes = [item for sublist in rows for item in sublist]

    singular_land = [np.array(ast.literal_eval(i), int) for i in list(set(landscapes))]

    path2 = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/all_data"
    files_data = os.listdir(path2)

    twos = set([nb2(i) for i in singular_land])
    nonzeros = set([nonzero(i) for i in singular_land])

    possible = list(itertools.product(nonzeros, twos))

    store = {}
    for i in range(len(possible)):
        store[possible[i]] = [j for j in singular_land if (nonzero(j) == possible[i][0] and nb2(j) == possible[i][1])]
    store2 = {k: v for k, v in store.items() if v}
    del store

    # region Matching & replacing procedure
    debut = time.time()
    replacement = {}

    for j in list(store2.keys()):
        if os.path.exists(path2 + "/land3_" + str(j[0]) + ".pkl"):
            filename = open(
                "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/all_data/land3_" + str(j[0]) + ".pkl", 'rb')
            set_to_check2 = pickle.load(filename)
            filename.close()
            if k == 0 and j == (0, 0):
                set_to_check = set_to_check2
            else:
                set_to_check = set()
                for subset in set_to_check2:
                    set_to_check.update(subset)

            for v in store2[j]:
                set_verif = utilities.equi_landscape(utilities.to_matrix(v))
                set_verif.add(tuple(v))
                list_verif = [*set_verif.intersection(set_to_check)]
                # new = list(list_verif[0])
                replacement[str(tuple(v))] = str(list_verif[0])

    for i in list(range(len(rows))):
        rows[i] = [replacement.get(n, n) for n in rows[i]]

    tester3 = pd.DataFrame(rows, dtype=str)
    tester3.to_csv(params.path + "/keys_succession/equiv_results_" + nums[2] + ".csv", index=False)

    return print('Dataset ' + str(params.files[k]) + " took " + str(time.time() - debut) + ' seconds')
