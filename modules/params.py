"""
Module for parameter set-up for Mouysset & Jean, 2023.

Parameters :
-----------
- Budget
- Landscape size (n) in nxn matrix
- Biodiversity values range for various cases
- Potential prescribed burns given budget.

Functions :
-------------------
- to_matrix returns matrix of landscape
- costs : returns homogeneous cost grid

"""

#region Packages
import os
import numpy as np
import itertools
#endregion

# Values
budget = 4
size = 3
R = size**2

# Biodiversity and fire thresholds
d_seuil = 2
d = d_seuil*np.ones(R)
m_seuil = 1
m = m_seuil*np.ones(R)
A_max = 2
nb_age = A_max+1
biodiv = np.arange(2, 101, 2)
biodiv3 = np.arange(0, 50, 2)
biodiv5 = np.arange(0,170,2)

# paths
if size == 4:
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_" + str(budget) + "/"

    files = os.listdir(path)
    files.remove('keys_succession')
    files.remove('successions_matched')
    files.remove('statistics')

elif size == 3:
    path_to_save_data = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/all_data"
    path_to_save_results = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results"

    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(budget) + "/"


    def structure(budget):
        sub_results = ['statistics', 'keys_succession', 'successions_matched', "successions_centrality",
                       "stat_centrality"]
        # Make directory to save data
        try:
            os.mkdir(path_to_save_data)
        except FileExistsError:
            pass
        # Make directory to save results
        try:
            os.mkdir(path_to_save_results)
        except FileExistsError:
            pass
        # Make directories inside each budget repository for specific intermediary and final products
        for b in range(1, budget):
            directory = "budget_" + str(b)
            path = os.path.join(path_to_save_results, directory)
            try:
                os.mkdir(path)
            except FileExistsError:
                pass

            for sub in sub_results:
                path2 = os.path.join(path, sub)
                try:
                    os.mkdir(path2)
                except FileExistsError:
                    pass

    structure(10)

    files = os.listdir(path)
    files.remove('keys_succession')
    files.remove('successions_matched')
    files.remove('statistics')
    files.remove('successions_centrality')
    files.remove('stat_centrality')


def to_matrix(land):
    """
    This function translates the array of length size into a sqrt(size)xsqrt(size) matrix
    :param land: np.array or tuple
    :return: matrix
    """
    if size == 3:
        return(np.array([list(land[0:3]),
                         list(land[3:6]),
                         list(land[6:9])]))
    elif size == 4:
        essai2 = np.array([list(land[0:4]),
                           list(land[4:8]),
                           list(land[8:12]),
                           list(land[12:16])])
        return essai2
    elif size == 5:
        essai2 = np.array([list(land[0:5]),
                           list(land[5:10]),
                           list(land[10:15]),
                           list(land[15:20]),
                           list(land[20:25])])
        return essai2
def cost():
    '''
    Generates the right cost matrix depending on size parameter
    :return: nd.array
    '''
    if size == 4:
        return(np.array([[1, 1, 1, 1],
                        [1, 1, 1, 1],
                        [1, 1, 1, 1],
                        [1, 1, 1, 1]]))
    elif size == 3:
        return(np.array([[1, 1, 1],
                         [1, 1, 1],
                         [1, 1, 1]]))
    elif size == 5:
        return (np.array([[1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1]]))

costs = cost()

# Generate the potential treatments, without budget consideration
pot_fire_value = (0, 1)
pot_fire = list(itertools.product(pot_fire_value, repeat=R))


# Adjust potential treatments for their admissibility : have to be affordable.
costs_param = "homogeneous"

if costs_param != 'homogeneous':
    pot_fire_budget = [x for x in pot_fire if sum(sum(to_matrix(x)*costs)) <= budget]
else:
    pot_fire_budget = [x for x in pot_fire if sum(x)<=budget]

#endregion