"""
Module to recover optimal succession of landscapes, Mousset & Jean, 2023

Functions :
----------
- equi_landscape_recover : recover which transformation has been applied to landscape
- recover : apply transformation to retrieve right landscape
- recover_full2 : chain application of recover
- succession_dataframe3 : recover optimal successions for 3x3 landscapes dataframes
- succession_dataframe4 : recover optimal successions for 4x4 landscapes dataframes

"""


import modules.utilities as utilities
import modules.params as params
import modules.matcher as matcher
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

def equi_landscape_recover(essai):
    """
    This function computes the unique matrices equivalent to the current 4x4 landscape using :
    - the initial matrix and 3 consecutive 90° rotations
    - A matrix transformed using horizontal symmetry and 3 consecutive 90° rotations
    - A matrix transformed using vertical symmetry
    """
    # 3 transformations
    m2 = np.rot90(essai)
    m3 = np.rot90(m2)
    m4 = np.rot90(m3)
    # Horizontal symmetry

    if params.size == 3:
        hor_sym = [essai[2], essai[1], essai[0]]
    elif params.size == 4:
        hor_sym = [essai[3], essai[2], essai[1], essai[0]]
    elif params.size == 5:
        hor_sym = [essai[4], essai[3], essai[2], essai[1], essai[0]]

    # 3 clockwise transformations
    m6 = np.rot90(hor_sym)
    m7 = np.rot90(m6)
    m8 = np.rot90(m7)
    evalu = [tuple(list(itertools.chain(*m2))),
             tuple(list(itertools.chain(*m3))),
             tuple(list(itertools.chain(*m4))),
             tuple(list(itertools.chain(*hor_sym))),
             tuple(list(itertools.chain(*m6))),
             tuple(list(itertools.chain(*m7))),
             tuple(list(itertools.chain(*m8))),
             tuple(list(itertools.chain(*essai)))]
    return evalu

def recover(land, transfo):
    land = utilities.to_matrix(land)
    if transfo == 0:
        m2 = np.rot90(land, axes=(1, 0))
        return tuple(list(itertools.chain(*m2)))
    if transfo == 1:
        m3 = np.rot90(np.rot90(land, axes=(1, 0)), axes=(1, 0))
        return tuple(list(itertools.chain(*m3)))
    if transfo == 2:
        m4 = np.rot90(np.rot90(np.rot90(land, axes=(1, 0)), axes=(1, 0)), axes=(1, 0))
        return tuple(list(itertools.chain(*m4)))
    if transfo == 3:
        if params.size == 4:
            m5 = [land[3], land[2], land[1], land[0]]
        elif params.size == 3:
            m5 = [land[2], land[1], land[0]]
        elif params.size == 5:
            m5 = [land[4], land[3], land[2], land[1], land[0]]
        return tuple(list(itertools.chain(*m5)))
    if transfo == 4:
        m6 = np.rot90(land, axes=(1, 0))
        if params.size == 4:
            m6 = [m6[3], m6[2], m6[1], m6[0]]
        elif params.size == 3:
            m6 = [m6[2], m6[1], m6[0]]
        elif params.size == 5:
            m6 = [m6[4], m6[3], m6[2], m6[1], m6[0]]
        return tuple(list(itertools.chain(*m6)))
    if transfo == 5:
        m7 = np.rot90(np.rot90(land, axes=(1, 0)), axes=(1, 0))
        if params.size == 4:
            m7 = [m7[3], m7[2], m7[1], m7[0]]
        elif params.size == 3:
            m7 = [m7[2], m7[1], m7[0]]
        elif params.size == 5:
            m7 = [m7[4], m7[3], m7[2], m7[1], m7[0]]
        return tuple(list(itertools.chain(*m7)))
    if transfo == 6:
        m8 = np.rot90(np.rot90(np.rot90(land, axes=(1, 0)), axes=(1, 0)), axes=(1, 0))
        if params.size == 4:
            m8 = [m8[3], m8[2], m8[1], m8[0]]
        elif params.size == 5:
            m8 = [m8[4], m8[3], m8[2], m8[1], m8[0]]
        elif params.size == 3:
            m8 = [m8[2], m8[1], m8[0]]
        return tuple(list(itertools.chain(*m8)))
    if transfo == 7:
        return tuple(list(itertools.chain(*land)))
    if transfo == "None":
        return tuple(list(itertools.chain(*land)))

def recover_full2(land, transfo, index):
    """
    Uses a list of transformations and the index of the land to recover the exact land succession from matched data.
    :param land: np.array
    :param transfo: int
    :param index: int
    :return: tuple
    """
    land3 = land
    if index == 0:
        return land3
    if index == 1:
        return recover(land3, transfo[0])
    if index == 2:
        return recover(recover(land3, transfo[1]), transfo[0])
    if index == 3:
        return recover(recover(recover(land3, transfo[2]), transfo[1]), transfo[0])
    if index == 4:
        return recover(recover(recover(recover(land3, transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 5:
        return recover(recover(recover(recover(recover(land3, transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 6:
        return recover(recover(recover(recover(recover(recover(land3, transfo[5]), transfo[4]), transfo[3]), transfo[2]),transfo[1]), transfo[0])
    if index == 7:
        return recover(recover(recover(recover(recover(recover(recover(land3, transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 8:
        return recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 9:
        return recover(recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[8]), transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 10:
        return recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[9]), transfo[8]), transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 11:
        return recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[10]), transfo[9]), transfo[8]), transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 12:
        return recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[11]), transfo[10]), transfo[9]), transfo[8]), transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 13:
        return recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[12]), transfo[11]), transfo[10]), transfo[9]), transfo[8]), transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 14:
        return recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[13]), transfo[12]), transfo[11]), transfo[10]), transfo[9]), transfo[8]), transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 15:
        return recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[14]), transfo[13]), transfo[12]), transfo[11]), transfo[10]), transfo[9]), transfo[8]), transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 16:
        return recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[15]), transfo[14]), transfo[13]), transfo[12]), transfo[11]), transfo[10]), transfo[9]), transfo[8]), transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 17:
        return recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[16]), transfo[15]), transfo[14]), transfo[13]), transfo[12]), transfo[11]), transfo[10]), transfo[9]), transfo[8]), transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 18:
        return recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[17]), transfo[16]), transfo[15]), transfo[14]), transfo[13]), transfo[12]), transfo[11]), transfo[10]), transfo[9]), transfo[8]), transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])
    if index == 19:
        return recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3, transfo[18]), transfo[17]), transfo[16]), transfo[15]), transfo[14]), transfo[13]), transfo[12]), transfo[11]), transfo[10]), transfo[9]), transfo[8]), transfo[7]), transfo[6]), transfo[5]), transfo[4]), transfo[3]), transfo[2]), transfo[1]), transfo[0])

def succession_dataframe3(data_path, biodiv_index, time_param='Yes', follow_up="No"):
    """
    Recover optimal successions for 20 periods and biodiversity constraint level for 3x3 landscapes
    :param data_path: str
        path to data
    :param biodiv_index: int
        biodiversity constraint level
    :param time_param: str
        optional, for time control
    :param follow_up: str
        optional, shows data header
    :return: .csv
        dataset of successions
    """
    debut = time.time()
    #region Load the current data

    file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(params.budget) + "/"+data_path)
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows_data = []
    for row in csvreader:
        row_add = [row[0], row[int(biodiv_index/2)], row[int(51+biodiv_index/2)]]
        rows_data.append(row_add)
    #endregion
    #region Import data to scan through - dict format and list format
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(params.budget) + "/"
    files = os.listdir(path)
    files.remove('keys_succession')
    files.remove('successions_matched')
    files.remove('statistics')
    files.remove('stat_centrality')
    files.remove('successions_centrality')
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[2]) in list(range(9 - params.budget, 10))]
    datas = {}
    datas[data_path] = rows_data
    for x in files:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(params.budget) + "/" + x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = [row[0], row[int(biodiv_index/2)], row[int(51+biodiv_index/2)]]
            rows_dict.append(row_add)
        datas[x] = rows_dict
    # endregion

    # region Import equivalent data - dict format and list formats
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(params.budget) + "/keys_succession/"
    files = os.listdir(path)
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[0]) in list(range(9 - params.budget, 10))]
    datas_equiv = {}
    for x in files:
        # region Load data in list format
        file = open(path + x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2)]
            rows_dict.append(row_add)
        datas_equiv[x] = rows_dict
        #endregion
    # endregion
    #region Add current data to dict repositories

    match_equiv = re.findall(r'\d+', str(data_path))
    file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(params.budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+".csv")
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows_dict = []
    for row in csvreader:
        row_add = row[int(biodiv_index/2)]
        rows_dict.append(row_add)
    datas_equiv[('equiv_results_'+str(match_equiv[2])+'.csv')] = rows_dict

    #endregion

    rows_out = []

    for p in list(range(len(rows_data))):
    # input
        land = str(rows_data[p][0])
    # output
        succession_list = [land]
        transformation = ['None']
        equiv_key = [land]
        value = [np.nan]
    #Meta parameters
        path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(params.budget) + "/"
        step = 0
        while step <= 20:
            check = matcher.nonzero(np.array(ast.literal_eval(equiv_key[-1])))
            #print(check, check2)
    # use relevant result data from loaded data :
            namer = ("land3_budget_" + str(params.budget) + "_" + str(check) + ".csv")

            if os.path.exists(path + namer):
                if namer in list(datas.keys()):
                    d2 = datas[namer]
                else:
                    file = open(path + "land3_budget_" + str(params.budget) + "_" + str(check) + ".csv")
                    csvreader = csv.reader(file)
                    header = next(csvreader)

                    rows_dict = []
                    for row in csvreader:
                        row_add = [row[0], row[int(biodiv_index/2)], row[int(51+biodiv_index/2)]]
                        rows_dict.append(row_add)
                    datas["land3_budget_" + str(params.budget) + "_" + str(check) + ".csv"] = rows_dict
                    d2 = rows_dict

                for g in list(range(len(d2))):
                    if d2[g][0] == equiv_key[-1]:
                        succession_index = g

                name2 = "equiv_results_" + str(check) + ".csv"

            succession_key = d2[succession_index][1]
            value_column = int(51 + biodiv_index / 2)

            if name2 in list(datas_equiv.keys()):
                d_equiv = datas_equiv[name2]
            else:
                file = open(path + "keys_succession/"+name2)
                csvreader = csv.reader(file)
                header = next(csvreader)

                rows_dict = []
                for row in csvreader:
                    row_add = row[int(biodiv_index/2)]
                    rows_dict.append(row_add)
                d_equiv = rows_dict
                datas_equiv[name2] = d_equiv

            #region Computations
            candidate_succ = d_equiv[succession_index]
            mat_succession = utilities.to_matrix(ast.literal_eval(succession_key))
            transfo = equi_landscape_recover(mat_succession).index(ast.literal_eval(candidate_succ))
            #endregion
            #region  list updates
            value.append(d2[succession_index][2])
            succession_list.append(succession_key)
            equiv_key.append(candidate_succ)
            transformation.append(transfo)
            #endregion

            step += 1
        #region Final processing : recovery of optimal succession with chain transformation
        c = list(map(ast.literal_eval, succession_list))
        final_succession = []
        for i in list(range(20)):
            e = transformation[0: i]
            # deleted +1
            final_succession.append(recover_full2(c[i], transfo=e, index=i))

        final_succession = final_succession + value
        rows_out.append(final_succession)
        if follow_up == 'Yes':
            print(str(100*len(rows_out) / len(rows_data))+" %")
        #endregion

    data_path_saved = data_path.replace('.csv', '_biodiv_index_' + str(biodiv_index) + '.csv')
    data_return = pd.DataFrame(rows_out)
    data_return.to_csv(path + 'successions_matched/' + data_path_saved)


    if time_param == 'Yes':
        return data_path + " is done and took " + str(time.time() - debut) + ' seconds'
    else :
        return data_path + " is done"

def succession_dataframe(data_path, biodiv_index,time_param='Yes',follow_up="No"):
    """
    Recover optimal successions for 20 periods and biodiversity constraint level for 3x3 landscapes
    :param data_path: str
        path to data
    :param biodiv_index: int
        biodiversity constraint level
    :param time_param: str
        optional, for time control
    :param follow_up: str
        optional, shows data header
    :return: .csv
        dataset of successions
    """
    debut = time.time()
    #region Load the current data

    file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_" + str(params.budget) + "/"+data_path)
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows_data = []
    for row in csvreader:
        row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
        rows_data.append(row_add)
    #endregion
    #region Import data to scan through - dict format and list format
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_" + str(params.budget) + "/"
    files = os.listdir(path)
    files.remove('keys_succession')
    files.remove('successions_matched')
    try:
        files.remove("list_update.csv")
    except ValueError:
        pass

    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[2]) in list(range(16 - params.budget, 17))]
    datas = {}
    datas[data_path]=rows_data
    for x in files:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_" + str(params.budget) + "/" +x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
            rows_dict.append(row_add)
        datas[x] = rows_dict
    # endregion

    # region Import equivalent data - dict format and list formats
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_" + str(params.budget) + "/keys_succession/"
    files = os.listdir(path)
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[0]) in list(range(16 - params.budget, 17))]
    datas_equiv = {}
    for x in files:
        # region Load data in list format
        file = open(path + x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[x] = rows_dict
        #endregion
    # endregion
    #region Add current data to dict repositories

    match_equiv = re.findall(r'\d+', str(data_path))
    if len(match_equiv)>4:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_" + str(params.budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+".csv")
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[('equiv_results_'+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+'.csv')]=rows_dict
    else:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_" + str(params.budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+".csv")
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[('equiv_results_'+str(match_equiv[2])+"_"+str(match_equiv[3])+'.csv')]=rows_dict

    #endregion

    rows_out = []

    for p in list(range(len(rows_data))):
    # input
        land = str(rows_data[p][0])
    # output
        succession_list = [land]
        transformation = ['None']
        equiv_key = [land]
        value = [np.nan]
    #Meta parameters
        path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_" + str(params.budget) + "/"
        step = 0
        while step <= 20:
            check = matcher.nonzero(np.array(ast.literal_eval(equiv_key[-1])))
            check2 = matcher.nb2(np.array(ast.literal_eval(equiv_key[-1])))
            #print(check, check2)
    # use relevant result data from loaded data :
            namer = ("land4_budget_" + str(params.budget) + "_" + str(check) + "_" + str(check2) + ".csv")
            namer2 = ("land4_budget_" + str(params.budget) + "_" + str(check) + "_" + str(check2) + "_cut_0.csv")

            if os.path.exists(path + namer):
                if namer in list(datas.keys()):
                    d2 = datas[namer]
                else:
                    file = open(path + "land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv")
                    csvreader = csv.reader(file)
                    header = next(csvreader)

                    rows_dict = []
                    for row in csvreader:
                        row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
                        rows_dict.append(row_add)
                    datas["land4_budget_" + str(params.budget) + "_" + str(check) + "_" + str(check2) + ".csv"] = rows_dict
                    d2 = rows_dict

                for g in list(range(len(d2))):
                    if d2[g][0]==equiv_key[-1]:
                        succession_index=g

                name2 = "equiv_results_" + str(check) + "_" + str(check2) + ".csv"

            elif os.path.exists(path + namer2):
                if namer2 in list(datas.keys()):
                    to_load = [x for x in list(datas.keys()) if x.startswith('land4_budget_' + str(params.budget) + "_" + str(check) + "_" + str(check2))]
                else :
                    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_" + str(
                        params.budget) + "/"
                    files = os.listdir(path)
                    files.remove('keys_succession')
                    files.remove('successions_matched')
                    to_load = [x for x in files if
                               x.startswith('land4_budget_' + str(params.budget) + "_" + str(check) + "_" + str(check2))]
                    for y in to_load:
                        file = open(path + y)
                        csvreader = csv.reader(file)
                        header = next(csvreader)

                        rows_dict = []
                        for row in csvreader:
                            row_add = [row[0], row[int(biodiv_index / 2)], row[int(82 + biodiv_index / 2)]]
                            rows_dict.append(row_add)
                        datas[y] = rows_dict

                for x in to_load:
                    d2 = datas[x]
                    for g in list(range(len(d2))):
                        if d2[g][0]==equiv_key[-1]:
                            succession_index = g
                            indexes = re.findall(r'\d+', str(x))
                            keep = x
                            name2 = "equiv_results_" + str(check) + "_" + str(check2) + "_cut_" + str(
                                indexes[4]) + ".csv"

                d2=datas[keep]
            succession_key = d2[succession_index][1]
            value_column = int(82 + biodiv_index / 2)

            if name2 in list(datas_equiv.keys()):
                d_equiv=datas_equiv[name2]
            else:
                file = open(path + "keys_succession/"+name2)
                csvreader = csv.reader(file)
                header = next(csvreader)

                rows_dict = []
                for row in csvreader:
                    row_add = row[int(biodiv_index/2-1)]
                    rows_dict.append(row_add)
                d_equiv = rows_dict
                datas_equiv[name2] = d_equiv

            #region Computations
            candidate_succ = d_equiv[succession_index]
            mat_succession = utilities.to_matrix(ast.literal_eval(succession_key))
            transfo = equi_landscape_recover(mat_succession).index(ast.literal_eval(candidate_succ))
            #endregion
            #region  list updates
            value.append(d2[succession_index][2])
            succession_list.append(succession_key)
            equiv_key.append(candidate_succ)
            transformation.append(transfo)
            #endregion

            step += 1
        #region Final processing : recovery of optimal succession with chain transformation
        c = list(map(ast.literal_eval, succession_list))
        final_succession = []
        for i in list(range(20)):
            e = transformation[0: i ]
            # deleted +1
            final_succession.append(recover_full2(c[i], transfo=e, index=i))

        final_succession = final_succession + value
        rows_out.append(final_succession)
        if follow_up == 'Yes':
            print(str(100*len(rows_out) / len(rows_data))+" %")
        #endregion

    data_path_saved = data_path.replace('.csv', '_biodiv_index_' + str(biodiv_index) + '.csv')
    data_return = pd.DataFrame(rows_out)
    data_return.to_csv(path+'successions_matched/'+data_path_saved)


    if time_param == 'Yes':
        return(data_path + " is done and took "+str(time.time() - debut)+' seconds' )
    else :
        return(data_path + " is done")