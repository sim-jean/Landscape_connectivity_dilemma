'''
Experiment for large scale landscapes with prioritization algorithm

'''

import numpy as np
import itertools
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import math
import modules.params as params

#region Functions
def connectivity_matrix_nx(size):
    G = nx.grid_2d_graph(size, size)
    for edge in G.edges:
        G.edges[edge]['weight'] = 1
    G.add_edges_from([
        ((x, y), (x+1, y+1))
        for x in range(size-1)
        for y in range(size-1)
    ] + [
        ((x+1, y), (x, y+1))
        for x in range(size-1)
        for y in range(size-1)
    ], weight=1)
    A = nx.to_numpy_array(G)
    np.fill_diagonal(A,1)
    A = A.astype(int)
    return A

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
    size = int(math.sqrt(len(landscape)))
    nodes = []

    #Population dynamics
        # Keep a follow up landscape for future reference.
    land_post = np.minimum(land + np.repeat(1,len(land)),2)
    land_follow = np.minimum(land + np.repeat(1,len(land)),2)

    # Loop on all cells
    for step in range(size**2):
        #region Set graph of fuel
            #Set of nodes
        land2 = (land_follow == 2)
        land2 = land2.astype(int)
            #Adjacency matrix
        adj = []
        for i in range(size**2):
            #
            adj.append(land2[i] * connectivity_matrix_nx(size)[i] * land2)
        adj = np.array(adj)
            #Graph using networkx
        G = nx.from_numpy_matrix(adj)
        #endregion

        #region Locate and assign treatment with eigenvector centrality (take most central value)
        cent = nx.eigenvector_centrality(G, max_iter=100000)
        candidate_cell = max(cent, key=cent.get)
        #endregion

        #region Biodiversity habitat : score, adjacency, loss
            # Compute land biodiversity : values between 1 and 2
        land_biod = ((land_follow > 0) & (land_follow < 3))
        land_biod = land_biod.astype(int)

            # Compute score of landscape
        score = land_biod.dot(connectivity_matrix_nx(size)).dot(np.transpose(land_biod))

            # Adjacency matrix for biodiversity graph
        adj_biod = []
        for i in range(size**2):
            adj_biod.append(land_biod[i]*connectivity_matrix_nx(size)[i]*land_biod)

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
        if budget ==0 or len(nodes) == size**2:
            break
            #endregion

    # Eventually, wrap up assignment using land_follow
    cells_to_burn = np.where(land_follow == 0)[0]
    if len(cells_to_burn)>0:
        land_post[cells_to_burn] = 0
    return land_post

#endregion

#region Experiment with large landscapes

#Define size of grid
size = 12

#Define randome landscape
landscape = np.random.randint(0,3,size**2)

# Define budget as a share of potential treatment, assuming homogeneous costs
budget = int(round(0.22 * size**2))
if budget <= 5:
    budget = 5

# reference point for biodiversity value:
land_biod = np.repeat(2,len(landscape))
biod_ref = land_biod.dot(connectivity_matrix_nx(size)).dot(land_biod)


# Initiate sequence over 10 periods
landscapes = [landscape]
for x in range(10):
    a = treatment_allocation(landscapes[-1],2, budg_spec=budget)
    #a_mat = a.reshape((size, size))
    landscapes.append(a)

for landscape in landscapes:
    ax = plt.axes()
    a_mat = landscape.reshape(size,size)
    sns.heatmap(a_mat, ax=ax, linewidth=0.5, cmap="Greens", vmin=min(landscape), vmax=2, center=1)
    plt.show()

#endregion