'''
Features the functions used to realize the dynamic programming part of Mouysset & Jean, 2023.

Contains :
- data import functions for 3x3 landscapes and 4x4 landscapes : data_list2, data_list
- Unit level dynamic programming for each landscape, in 3x3 and 4x4 : lowlev3_dynprog, lowlev_dynprog
'''

import modules.utilities as utilities
import modules.params as params
import modules.outputs as outputs
#region Packages

import numpy as np
import os
import pickle
import pandas as pd
import networkx as nx
import time

#endregion

#region Data import
def data_list2(input):
    """
    Utility function for data import from 4x4 landscapes

    :param input: int
    :return: list
    """
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data/"
    files = os.listdir(path)

    filename = open(path+files[input], 'rb')
    list_nz = pickle.load(filename)
    filename.close()
    list_nz = list(list_nz)
    return list_nz

def data_list3(input):
    """
    Utility function for data import from 3x3 landscapes
    :param input: int
    :return: list
    """
    path = params.path_to_save_data
    files = os.listdir(path)

    filename = open(path +"/" +files[input], 'rb')
    list_nz = pickle.load(filename)
    filename.close()

    list_nz = [item for sublist in list_nz for item in sublist]

    return list_nz
#endregion


def lowlev3_dynprog(land):
    """
    Function that computes the value function from Mouysset & Jean, 2023 on 3x3 landscapes, and returns the optimal landscape successions, the min
    value and the optimal prescribed burns

    :param land: np.array
    :return: list - lands, prescribed burns, values
    """
    # Load the environment
    #data_list3(input)
    land = tuple(land)
    # Set up the output variables

    # compute the new land and associated value
    new_land = utilities.fuel_dyn(land, params.pot_fire_budget)
    value = utilities.high_fuel_con_vect(new_land) + utilities.high_fuel_connect(land) + 1000000 * (utilities.high_biod_con_vect(new_land) <= params.biodiv3)

    # find min
    a = value.min(0)
    b = value.argmin(0)

    store = [tuple(x) for x in new_land[np.asarray(b, int)][0].tolist()]
    listed_values= [[land], store, np.asarray(b)[0].tolist(), np.asarray(a)[0].tolist()]
    list_values = [item for sublist in listed_values for item in sublist]

    return list_values

def lowlev_dynprog(land):
    """
    Function that computes the value function from Mouysset & Jean, 2023 on 4x4 landscapes, and returns the optimal landscape successions, the min
    value and the optimal prescribed burns

    :param land: np.array
    :return: list - lands, prescribed burns, values
    """
    # Load the environment
    #data_list2(input)
    land = tuple(land)
    # Set up the output variables
    # compute the new land and associated value
    new_land = utilities.fuel_dyn(land, params.pot_fire_budget)
    value = utilities.high_fuel_con_vect(new_land) + utilities.high_fuel_connect(land) + 1000000 * (utilities.high_biod_con_vect(new_land) <= params.biodiv)
    # find min
    a = value.min(0)
    b = value.argmin(0)

    store = [tuple(x) for x in new_land[np.asarray(b, int)][0].tolist()]
    listed_values = [[land], store, np.asarray(b)[0].tolist(), np.asarray(a)[0].tolist()]
    list_values = [item for sublist in listed_values for item in sublist]

    return list_values



def treatment_allocation(land,biod_constraint,budg_spec=-1):
    """
    Algorithm for treatment allocation using eigencentrality
    :param land: np.array
        landscape
    :param biod_constraint: int
        Level of biodiv constraint
    :param budg_spec: int
        Potential for new budget
    :return:
    """
    if budg_spec>0:
        budget = budg_spec
    else:
        budget = params.budget

    nodes = []

    #Population dynamics
        # Keep a follow up landscape for future reference.
    land_post = np.minimum(land + np.repeat(1,len(land)),2)
    land_follow = np.minimum(land + np.repeat(1,len(land)),2)

    # Loop on all cells
    for step in range(params.size**2):
        #region Set graph of fuel
            #Set of nodes
        land2 = (land_follow == 2)
        land2 = land2.astype(int)
            #Adjacency matrix
        adj = []
        for i in range(params.size**2):
            #
            adj.append(land2[i] * utilities.connectivity_mat()[i] * land2)
        adj = np.array(adj)
            #Graph using networkx
        G = nx.from_numpy_matrix(adj)
        #endregion

        #region Locate and assign treatment with eigenvector centrality (take most central value)
        cent = nx.eigenvector_centrality(G, max_iter=1000)
        candidate_cell = max(cent, key=cent.get)
        #endregion

        #region Biodiversity habitat : score, adjacency, loss
            # Compute land biodiversity : values between 1 and 2
        land_biod = ((land_follow > 0) & (land_follow < 3))
        land_biod = land_biod.astype(int)

            # Compute score of landscape
        score = land_biod.dot(utilities.connectivity_mat()).dot(np.transpose(land_biod))

            # Adjacency matrix for biodiversity graph
        adj_biod = []
        for i in range(params.size**2):
            adj_biod.append(land_biod[i]*utilities.connectivity_mat()[i]*land_biod)

        adj_biod = np.array(adj_biod)
            # Compute biodiversity loss :
        biod_loss = 2 * (sum(adj_biod[candidate_cell])-1) + 1
        #endregion

        #region Treatment assignment or alternative
        if score - biod_loss >= biod_constraint and budget >0:
            # If biodiversity remains sufficient, assign treatment and update budget

            land_follow[candidate_cell] = 0
            budget -= 1
            nodes.append(0)
            # If biodiversity score is not satisfied, rule out this cell and start again.
        else :
            land_follow[candidate_cell] = 10
            nodes.append(0)
            # Other conditions : either budget is 0, or all nodes have been scanned and none was suited for treatment assignement
        if budget ==0 or len(nodes) == params.size**2:
            break
            #endregion

    # Eventually, wrap up assignment using land_follow
    cells_to_burn = np.where(land_follow == 0)[0]
    if len(cells_to_burn)>0:
        land_post[cells_to_burn] = 0
    return land_post

def land_succ(lands):
    """
    Generates 20 period succession for initial landscapes for all biodiversity levels using treatment assignment based on eigenvector centrality
    :param lands: np.array
        landscapes
    :return: data, .csv
    """
    # Initiate storage
    storage = {}
    storage['index'] = []
    storage['biod'] = []
    storage['value'] = []

    start = time.time()
    index = 0

    # Pick biodiversity range depending on landscape size
    if params.size==3:
        biodiv_range = list(range(0,max(params.biodiv3),2))
    elif params.size==4:
        biodiv_range = list(range(0, max(params.biodiv), 5))

    # Generate the vector of biodiversity to add to data each step
    biodiv_coord = [[x]*20 for x in biodiv_range]
    biodiv_coord = [item for sublist in biodiv_coord for item in sublist]

    # Loop over land
    for land in lands :
        succ = []
        # Loop over biodiversity constraints
        for biod in biodiv_range:
            # Initial land and successions resulting from treatment allocation
            successions = [land]
            for t in range(19):
                successions.append(treatment_allocation(successions[-1],biod))
            succ.extend(successions)

        # Store data
        storage['index'].extend([index]*len(succ))
        storage['biod'].extend(biodiv_coord)
        storage['value'].extend(succ)
        index += 1
        #print(index/len(lands))
    print("File took :" + str(time.time()-start))

    # Add budget and save
    storage['budget']=[params.budget]*len(storage['value'])
    check = pd.DataFrame(storage)
    return check

