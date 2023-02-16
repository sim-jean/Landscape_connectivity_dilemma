"""
Module with utility functions for Mouysset & Jean, 2O23

Functions :
-----------
- structure : creates folder structure up to maximum budget specified

- to_matrix : returns matrix form for landscape
- equi_landscape : returns set of equivalent landscapes
- equi_landscape_recover : returns set of equivalent landscapes and self for recovery

- low_lev_land3 : generate 3x3 landscapes
- low_lav_land4 : generate 4x4 landscapes

- connectivity_mat : connectivity matrix for different sizes
- fuel_dyn : growth function for landscape
- high_fuel : indicator for high fuel or not in a landscape
- high_fuel_connect : score for fuel connectivity
- high_fuel_con_vect : vectorized version of it
- high_biod : indicator for biodiversity habitat
- high_biod_con_vect : score for biodiversity habitat connectivity, vectorized
"""

#region Packages
import modules.params as params
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


#region Folders
#endregion
#region Equivalent landscapes : self excluded and included
def to_matrix(land):
    """
    This function translates the array of length size into a sqrt(size)xsqrt(size) matrix
    :param land: np.array or tuple
    :return: matrix
    """
    if params.size == 3:
        return(np.array([list(land[0:3]),
                         list(land[3:6]),
                         list(land[6:9])]))
    elif params.size == 4:
        essai2 = np.array([list(land[0:4]),
                           list(land[4:8]),
                           list(land[8:12]),
                           list(land[12:16])])
        return essai2
    elif params.size == 5:
        essai2 = np.array([list(land[0:5]),
                           list(land[5:10]),
                           list(land[10:15]),
                           list(land[15:20]),
                           list(land[20:25])])
        return essai2
def equi_landscape(land):
    """
    This function computes the unique matrices equivalent to the current 4x4 landscape using :
    - the initial matrix and 3 consecutive 90° rotations
    - A matrix transformed using horizontal symmetry and 3 consecutive 90° rotations
    - A matrix transformed using vertical symmetry
    :param land: matrix
    :return: list of tuples
    """
    # 3 transformations clockwise rotation
    m2 = np.rot90(land)
    m3 = np.rot90(m2)
    m4 = np.rot90(m3)
    # Horizontal symmetry
    if params.size == 4:
        hor_sym = [land[3], land[2], land[1], land[0]]
    elif params.size == 3:
        hor_sym = [land[2], land[1], land[0]]
    elif params.size==5:
        hor_sym = [land[4], land[3], land[2], land[1], land[0]]
    # 3 clockwise transformations
    m6 = np.rot90(hor_sym)
    m7 = np.rot90(m6)
    m8 = np.rot90(m7)
    evalu = (tuple(list(itertools.chain(*m2))),
             tuple(list(itertools.chain(*m3))),
             tuple(list(itertools.chain(*m4))),
             tuple(list(itertools.chain(*hor_sym))),
             tuple(list(itertools.chain(*m6))),
             tuple(list(itertools.chain(*m7))),
             tuple(list(itertools.chain(*m8))))
    return set(evalu)

def equi_landscape_recover(land):
    """
    Returns the same results as equi_landscape, with self added
    :param land: matrix
    :return: list of tuples
    """
    # Rotate clockwise
    m2 = np.rot90(land)
    m3 = np.rot90(m2)
    m4 = np.rot90(m3)

    # Define horizontal symmetry
    if params.size == 3:
        hor_sym = [land[2], land[1], land[1]]
    elif params.size == 4:
        hor_sym = [land[3], land[2], land[1], land[0]]
    elif params.size == 5:
        hor_sym = [land[4], land[3], land[2], land[1], land[0]]

    # Rotate clockwise
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
             tuple(list(itertools.chain(*land)))]
    return evalu
#endregionn
#region Generate unique landscapes
def low_lev_land3(nb2, non_z):
    """
    Returns the unique landscapes from a combination of non zero elements and a number of twos. Uses a nested combinatoric analysis : first, 0 and non zero.
    Then, among non zeros, elements take value 1 and 2
    :param nb2: int
        Number of desired twos
    :param non_z: int
        Number of desired non zeros
    :return: set
        Set of unique landscapes of size 3x3
    """
    unique_landscapes_nb2 = set()
    size = 9
    non_zero = [1, 2]
    # Non zero elements are 1 or 2

    # Following three lines generate potential combinations of where non zero elements are located, in the form of a 0-1 vector.
    which = np.array(list(itertools.combinations(range(size), non_z)))
    grid = np.zeros((len(which), size), dtype="int8")
    grid[np.arange(len(which))[None].T, which] = 1

    # Loop on each potential arrangement of 0 and non-zero
    for i in np.array(range(len(grid))):
        grid_arr = np.array(grid[i])

        # Find non zero elements to be replaced
        index_replace = np.where(grid_arr == 1)

        # Generate all the potential combinations of 1 and 2 of length nonzero
        candidat_lists = list(itertools.product(non_zero, repeat=non_z))

        # Restrict this list to candidates such that the number of 1 and 2 is specified by inputs
        nb2_list = [x for x in candidat_lists if sum(x) == (non_z-nb2)+nb2*2]

        # Loop on potential 1 and 2 arrangements among non-zeros
        for j in nb2_list:
            candidate = grid_arr

            # Replace 1 in grid_arr with combination of 2 and 1
            for k in np.array(range(len(j))):
                candidate[index_replace[0][k]] = j[k]

            # Generate all equivalent landscapes to this candidate
            mat_test = equi_landscape(to_matrix(candidate))

            # Update landscape list only if no equivalent is present in the list
            if len(mat_test.intersection(unique_landscapes_nb2)) == 0:
                unique_landscapes_nb2.add(tuple(candidate))
    return unique_landscapes_nb2

def low_lev_land4(nb2, non_z):
    """
    Returns the unique landscapes from a combination of non zero elements and a number of twos. Uses a nested combinatoric analysis : first, 0 and non zero.
    Then, among non zeros, elements take value 1 and 2
    :param nb2: int
        Number of desired twos
    :param non_z: int
        Number of desired non zeros
    :return: set
        Set of unique landscapes of size 4x4

    """
    unique_landscapes_nb2 = set()
    size = 16
    non_zero = [1, 2]
    # Non zero elements are 1 or 2

    # Following three lines generate potential combinations of where non zero elements are located, in the form of a 0-1 vector.
    which = np.array(list(itertools.combinations(range(size), non_z)))
    grid = np.zeros((len(which), size), dtype="int8")
    grid[np.arange(len(which))[None].T, which] = 1

    # Loop on each potential arrangement of 0 and non-zero
    for i in np.array(range(len(grid))):
        grid_arr = np.array(grid[i])

        # Find non zero elements to be replaced
        index_replace = np.where(grid_arr == 1)

        # Generate all the potential combinations of 1 and 2 of length nonzero
        candidat_lists = list(itertools.product(non_zero, repeat=non_z))

        # Restrict this list to candidates such that the number of 1 and 2 is specified by inputs
        nb2_list = [x for x in candidat_lists if sum(x) == (non_z-nb2)+nb2*2]

        # Loop on potential 1 and 2 arrangements among non-zeros
        for j in nb2_list:
            candidate = grid_arr

            # Replace 1 in grid_arr with combination of 2 and 1
            for k in np.array(range(len(j))):
                candidate[index_replace[0][k]] = j[k]

            # Generate all equivalent landscapes to this candidate
            mat_test = equi_landscape(to_matrix(candidate))

            # Update landscape list only if no equivalent is present in the list
            if len(mat_test.intersection(unique_landscapes_nb2)) == 0:
                unique_landscapes_nb2.add(tuple(candidate))
    return unique_landscapes_nb2
#endregion
#region Dynamics, biodiversity and fuel
def connectivity_mat():
    """
    Returns the adjacency matrix for a King's Graph of dimensions size**2
    :return: np.array matrix
    """
    if params.size == 3:
        return np.array([[1, 1, 0, 1, 1, 0, 0, 0, 0],
                          [1, 1, 1, 1, 1, 1, 0, 0, 0],
                          [0, 1, 1, 0, 1, 1, 0, 0, 0],
                          [1, 1, 0, 1, 1, 0, 1, 1, 0],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [0, 1, 1, 0, 1, 1, 0, 1, 1],
                          [0, 0, 0, 1, 1, 0, 1, 1, 0],
                          [0, 0, 0, 1, 1, 1, 1, 1, 1],
                          [0, 0, 0, 0, 1, 1, 0, 1, 1]])

    elif params.size == 4:
        return np.array([[1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                          [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
                          [1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0],
                          [0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0],
                          [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0],
                          [0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0],
                          [0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1],
                          [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1],
                          [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1]])

    elif params.size == 5:
        return np.array([[1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 1
                          [1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 2
                          [0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 3
                          [0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 4
                          [0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 5
                          [1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 6
                          [1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 7
                          [0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 8
                          [0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 9
                          [0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 10
                          [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # 11
                          [0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],  # 12
                          [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0],  # 13
                          [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0],  # 14
                          [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],  # 15
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0],  # 16
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0],  # 17
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],  # 18
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1],  # 19
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1],  # 20
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0],  # 21
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0],  # 22
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],  # 23
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1],  # 24
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1]])

    else:
        print("No connectivity matrix associated with this size")

def fuel_dyn(landscape, presc_burn):
    """
    Returns the fuel dynamics resulting from (i) natural growth and (ii) prescribed burns
    :param landscape: np.array or tuple
    :param presc_burn: np.array
    :return:
    """
    post = (np.array(landscape) + 1)*(1-np.array(presc_burn))
    land_post = np.minimum(post, 2)
    return land_post

def high_fuel(landscape):
    """
    This function returns a Boolean vector of fuel content larger than fuel risk threshold for each cell.
    :param landscape: np.array or tuple
    :return: boolean
    """
    return (np.array(landscape) >= params.d)

def high_fuel_connect(landscape):
    """
    Computes the fuel connectivity score for a given landscape.
    :param landscape: np.array or tuple
    :return: int
    """
    return(high_fuel(landscape).dot(connectivity_mat()).dot(high_fuel(landscape)[:, np.newaxis]))

def high_fuel_con_vect(land):
    """
    Vectorized version of high_fuel_connect
    :param land: np.array or tuple
    :return: int
    """
    return(np.transpose(np.matrix(np.diagonal((land >= params.d).dot(connectivity_mat()).dot(np.transpose((land >= params.d)))))))

def high_biodiv(landscape):
    """
    This function returns a Boolean vector of fuel content larger than the biodiversity value threshold for each cell.
    :param landscape: np.array or tuple
    :return: bool
    """
    return (np.array(landscape) >= params.m_seuil)

def high_biod_con_vect(new_land):
    """
    Equivalent of fuel_con_vect for biodiversity
    :param new_land: np.array or tuple
    :return: int
    """
    a = np.mat(np.diagonal((new_land >= params.m).dot(connectivity_mat()).dot(np.transpose((new_land >= params.m))))).T
    return a
#endregion