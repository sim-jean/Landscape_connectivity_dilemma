"""
Module for the outputs in Mouysset & Jean, 2023

Classes :
---------
- Class : Graph, Q_Item for graph analysis
- Class : Land for landscape analysis

Functions :
------------
- All the functions defined for classes
- convergence : for convergence of landscape series
- statistics_dataframe : returns statistics of interest computed for whole dataframes
- view3D : view the envelope of the optimization program.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import statistics as stats
import ast
import re
import pandas as pd
import time
import matplotlib
import scipy
import networkx as nx

import current_scripts.utilities as utilities
import current_scripts.params as params


class Graph:
    """
    Class used to represent graphs

    ...
    Attributes
    ----------
    row : int
        Row considered
    col : int
        Column considered
    g : list
        Graph associated

    Methods
    ----------
    isSafe()
        Verifies if can cell can be visited for DFS (in the matrix, and non visited)
    DFS()
        Performs DFS on matrix
    countIslands()
        Counts the number of non-adjacent components

    """
    def __init__(self, row, col, g):
        self.ROW = row
        self.COL = col
        self.graph = g

    def isSafe(self, i, j, visited):
        '''
        A function to check if a given cell (row, col) can be included in DFS
        :param i: int
        :param j: int
        :param visited: bool
        :return:
        '''
        # row number is in range, column number
        # is in range and value is 1
        # and not yet visited
        return (i >= 0 and i < self.ROW and
                j >= 0 and j < self.COL and
                not visited[i][j] and self.graph[i][j])

    def DFS(self, i, j, visited):
        '''
        A utility function to do DFS for a 2D boolean matrix. It only considers the 8 neighbours as adjacent vertices
        :param i: int
        :param j: int
        :param visited:  bool
        :return:
        '''

        # These arrays are used to get row and
        # column numbers of 8 neighbours
        # of a given cell
        rowNbr = [-1, -1, -1, 0, 0, 1, 1, 1]
        colNbr = [-1, 0, 1, -1, 1, -1, 0, 1]

        # Mark this cell as visited
        visited[i][j] = True

        # Recur for all connected neighbours
        for k in range(8):
            if self.isSafe(i + rowNbr[k], j + colNbr[k], visited):
                self.DFS(i + rowNbr[k], j + colNbr[k], visited)

    def countIslands(self):
        '''
        Returns the count of islands in graph from adjacency matrix
        :return: int
        '''
        # Make a bool array to mark visited cells.
        # Initially all cells are unvisited
        visited = [[False for j in range(self.COL)]for i in range(self.ROW)]

        # Initialize count as 0 and traverse
        # through the all cells of
        # given matrix
        count = 0
        area = [0]
        zone_l=[]
        for i in range(self.ROW):
            for j in range(self.COL):
                # If a cell with value 1 is not visited yet,
                # then new island found
                if visited[i][j] == False and self.graph[i][j] == 1:
                    # Visit all cells in this island
                    # and increment island count

                    self.DFS(i, j, visited)
                    count += 1
                    area.append(sum([sum(x) for x in visited])-area[-1])
                    zone_l.append((i, j))
        area.pop(0)

        return count, area, zone_l
class QItem:
    """
    A class used to represent points in matrix associated to a graph

    Attributes
    ----------
    row : int
        integer describing row
    col : int
        integer describing col
    distance : int
        integer describing some distance from point


    """
    def __init__(self, row, col, dist):
        self.row = row
        self.col = col
        self.dist = dist

    def __repr__(self):
        return f"QItem({self.row}, {self.col}, {self.dist})"


def minDistance(grid, sourcex, sourcey, destx, desty):
    """
    returns the minimum distance between two cells in an adjacency matrix
    :param grid: matrix
    :param sourcex: int
    :param sourcey: int
    :param destx: int
    :param desty: int
    :return: float
    """
    source = QItem(sourcex, sourcey, 0)


    # To maintain location visit status
    visited = [[False for _ in range(len(grid[0]))]
               for _ in range(len(grid))]

    # applying BFS on matrix cells starting from source
    queue = []
    queue.append(source)
    visited[source.row][source.col] = True
    while len(queue) != 0:
        source = queue.pop(0)

        # Destination found;
        if source.row == destx and source.col == desty:
            return source.dist

        # moving up left
        if isValid(source.row - 1, source.col - 1, grid, visited):
            queue.append(QItem(source.row - 1, source.col - 1, source.dist + math.sqrt(2)))
            visited[source.row - 1][source.col - 1] = True

        # moving up right
        if isValid(source.row - 1, source.col + 1, grid, visited):
            queue.append(QItem(source.row - 1, source.col + 1, source.dist + math.sqrt(2)))
            visited[source.row - 1][source.col + 1] = True

        # moving down right
        if isValid(source.row + 1, source.col + 1, grid, visited):
            queue.append(QItem(source.row + 1, source.col + 1, source.dist + math.sqrt(2)))
            visited[source.row + 1][source.col + 1] = True
        # moving down left
        if isValid(source.row + 1, source.col - 1, grid, visited):
            queue.append(QItem(source.row - 1, source.col + 1, source.dist + math.sqrt(2)))
            visited[source.row + 1][source.col - 1] = True

        # moving up
        if isValid(source.row - 1, source.col, grid, visited):
            queue.append(QItem(source.row - 1, source.col, source.dist + 1))
            visited[source.row - 1][source.col] = True

        # moving down
        if isValid(source.row + 1, source.col, grid, visited):
            queue.append(QItem(source.row + 1, source.col, source.dist + 1))
            visited[source.row + 1][source.col] = True

        # moving left
        if isValid(source.row, source.col - 1, grid, visited):
            queue.append(QItem(source.row, source.col - 1, source.dist + 1))
            visited[source.row][source.col - 1] = True

        # moving right
        if isValid(source.row, source.col + 1, grid, visited):
            queue.append(QItem(source.row, source.col + 1, source.dist + 1))
            visited[source.row][source.col + 1] = True

    return 1000
    #dealt with 0 if no path exists, but it is problematic ; if all patches have the same size
    # and no
# checking where move is valid or not
def isValid(x, y, grid, visited):
    '''
    Verifies if cell can be visited by DFS algorithm
    :param x: int
        row
    :param y: int
        col
    :param grid: mat
        associated adjacency matrix to graph
    :param visited: bool
    :return: bool
    '''
    if ((x >= 0 and y >= 0) and
            (x < len(grid) and y < len(grid[0])) and
            (grid[x][y] != 0) and (visited[x][y] == False)):
        return True
    return False

class Land:
    """
    A class used to represent landscapes

    Attributes
    -----------
    shape: np.array
        vector form of landscape
    size: int
        number of rows and columns in a square matrix representation
    mat : nd.array
        matrix representation of landscape
    biod : int
        Biodiversity habitat score computed from landscape and maximum connectivity matrix
    fuel : int
        Fuel score computed from landscape and maximum connectivity matrix
    biod_relative : float
        Relative score compared to maximum biodiversity connectivity
    fuel_relative : float
        Relative score compared to maximum fuel connectivity
    node_biod : int
        Number of nodes in the biodiversity habitat graph associated with landscape, also number of 1s & 2s in landscape
    node_fuel : int
        Number of nodes in the fuel graph associated with landscape, also number of 2s in landscape
    zero : int
        Number of empty patches in landscape
    habitat_area_ratio : float
        Ratio of effective to potential habitat in landscape
    fuel_area_ratio : float
        Ratio of effective to potential fuel in landscape

    Methods
    ---------
    view()
        Returns heatmap of the matrix form of the landscape
    view_succession()
        Returns heatmap of the optimal succession of the landscape considering biodiversity constraint = index
    equivalent()
        Returns the set of equivalent landscapes resulting from rotations and symmetry.
    fuel_dynamics()
        Returns the potential landscapes resulting from fuel dynamics considering potential prescribed burns
    succession()
        Returns the optimal succession of landscape considering various biodiversity constraints, and the corresponding value
    test_algo()
        Returns the dictionary of paths in a graph where nodes are cells with values larger or equal/equal to k
    shortest_path(var="biod")
        Returns the shortest path between two cells in adjacency matrix associated to graph.
        Default graph considered is biodiversity habitat. Alternative : 'fuel' or 'zeros'.
    IIC(var="biod", nc=True)
        Integral Index of Connectivity (Pascual Horta, Suara, 2006) associated to graph
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
        Default nc is True : give a finite value to disconnected cells, if False, discarded
    CPL(var="biod", nc=True)
        Characteristic path length associated to graph
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
        Default nc is True : give a finite value to disconnected cells, if False, discarded
    diameter(var="biod")
        Diameter (longest path) associated to graph
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    landscape_shape_index(var="biod")
        Landscape shape index for raster data, see Fragstat documentation (p. 115)
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    diagonal(direction = "NO-SE", var="biod")
        Length of shortest diagonal path
        Default direction is North-West to South-East ('NO-SE'), alternative ('SO-NE')
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    number_empty_neighbor(var="biod")
        Number of empty neighbors in 4 directions
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    number_zero_adjacent(var="biod")
        Number of empty neighbors in 8 directions
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    components(var="biod")
        Returns the number of components/subgraphs in landscape/graph associated with variable
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    components_area(var="biod")
        Returns the area of components/ numbers of nodes of subgraphs associated to specific variable
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    components_area_max(var="biod")
        Returns area of largest component/ number of nodes of largest subgraph associated with variable
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    components_perimeter(var="biod")
        Returns perimeter of components in landscape associated with variable
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    components_shape_index(var="biod")
        Analogous to Landscape Shape Index at the largest component scale
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    node_degree(var="biod")
        Returns degree corresponding to node
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    connectivity_correlation(var="biod", option_plot="No")
        Returns Pearson coefficient of correlation between node degree and neighbors' degrees.
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    perimeter(var="biod")
        Returns perimeter of overall landscape, depending on variable
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"
    perimeter_to_area(var="biod")
        Returns ratio of perimeter to area on the whole landscape, depending on variable
        Default var is "biod": graph considered is associated with biodiversity habitat. Alternative : "fuel"

    """
    def __init__(self, shape):
        self.shape = shape
        self.size = int(math.sqrt(len(shape)))
        self.mat = utilities.to_matrix(shape)
        self.biod = np.isin(shape, [-1, 1, 2]).dot(utilities.connectivity_mat()).dot(np.transpose(np.isin(shape, [-1, 1, 2])))
        self.biod_relative = self.biod/(np.repeat(1, len(self.shape)).dot(utilities.connectivity_mat()).dot(np.transpose(np.repeat(1, len(self.shape)))))
        self.fuel = (shape >= 2).dot(utilities.connectivity_mat()).dot(np.transpose((shape >= 2)))
        self.fuel_relative = self.fuel/(np.repeat(1, len(self.shape)).dot(utilities.connectivity_mat()).dot(np.transpose(np.repeat(1, len(self.shape)))))
        self.node_biod = sum(shape > 0)
        self.node_fuel = sum(shape == 2)
        self.zero = sum(shape == 0)
        self.habitat_area_ratio = sum(self.shape >= 1)/(self.size**2)
        self.fuel_area_ratio = sum(self.shape == 2)/(self.size**2)

    def view(self):
        """
        Returns heatmap of the matrix form of the landscape

        :return: sns plot
        """
        ax = plt.axes()
        ax.set_title('Initial biodiv = '+str(self.biod))
        sns.heatmap(self.mat, ax= ax, linewidth=0.5, cmap="Greens", vmin=min(self.shape), vmax=2, center=1)
        plt.show()

    def view_succession(self, index):
        """
        Returns heatmap of the optimal succession of the landscape considering biodiversity constraint = index

        :param index: int - Biodiversity constraint level
        :return: sns plot
        """
        if len(index) == 1:
            ax = plt.axes()
            ax.set_title('Biodiv constraint = ' + str(index[0]) + " with budget = " + str(params.budget))
            y = utilities.to_matrix(self.succession_land[index[0]])
            sns.heatmap(y, ax=ax, linewidth=0.5, cmap="Greens", vmin=0, vmax=2, center=1)
            plt.show()
        else:
            number_x = 0
            rows = len(index)//2
            for i in range(rows):
                for j in range(2):
                    ax = plt.subplot2grid((rows, 2), (i, j))
                    ax.title.set_text("Biodiv constraint =" + str(index[number_x]) + ' with budget =' + str(params.budget))
                    y = utilities.to_matrix(self.succession_land[index[number_x]])
                    ax = sns.heatmap(y, linewidth=0.5, cmap="Greens", vmin=0, vmax=2, center=1)
                    number_x += 1
            plt.show()

    def equivalent(self):
        """
        Returns the set of equivalent landscapes resulting from rotations and symmetry.

        :return: set of tuples
        """
        m2 = np.rot90(self.mat)
        m3 = np.rot90(m2)
        m4 = np.rot90(m3)
        if self.size == 5:
             # Horizontal symmetry
            hor_sym = [self.mat[4], self.mat[3], self.mat[2], self.mat[1], self.mat[0]]
        elif self.size == 4:
            hor_sym = [self.mat[3], self.mat[2], self.mat[1], self.mat[0]]
        elif self.size == 3:
            hor_sym = [self.mat[2], self.mat[1], self.mat[0]]
            # 3 clockwise transformations
        else:
            ValueError("Size is not supported")
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

    def fuel_dynamics(self, presc_burn):
        """
        Returns the potential landscapes resulting from fuel dynamics considering potential prescribed burns

        :param presc_burn: np.array or tuple
        :return: np.array
        """
        post = (np.array(self.shape) + 1) * (1 - np.array(presc_burn))
        land_post = np.minimum(post, 2)
        self.fuel_dynamics = land_post

    def succession(self, constraint_biodiv):
        """
        Returns the optimal succession of landscape considering various biodiversity constraints, and the corresponding value

        :param constraint_biodiv: int
            Belongs to [0,100]
        :return: np.array/tuple and float
        """
        if len(self.fuel_dynamics) > 1:
            biod = np.mat(np.diagonal((self.fuel_dynamics >= 1).dot(utilities.connectivity_mat()).dot(np.transpose((self.fuel_dynamics >= 1))) - (self.fuel_dynamics >= 1).sum(
                    axis=1))).T
            fuel = np.transpose(np.matrix(np.diagonal(
                (self.fuel_dynamics >= 2).dot(utilities.connectivity_mat()).dot(np.transpose((self.fuel_dynamics >= 2))) - (self.fuel_dynamics >= 2).sum(axis=1))))

            value = fuel + self.fuel + 1000*(biod <= constraint_biodiv)
            a = value.min(0)
            b = value.argmin(0)

            self.succession_land = {k: v for (k, v) in zip(constraint_biodiv, [np.array(x) for x in self.fuel_dynamics[np.asarray(b, int)][0].tolist()])}
            self.succession_value ={k: v for (k, v) in zip(constraint_biodiv, np.asarray(a)[0].tolist())}
            return self.succession_land, self.succession_value

        else:
            print("Problem")
    # Functions for characterization

    def test_algo(self, k):
        '''
        Returns the dictionary of paths in a graph where nodes are cells with values larger or equal/equal to k

        :param k: int - belongs to {0,1,2}
        :return: dictionary
        '''
        if k >= 0:
            grid = (self.mat >= k)
        elif k == 0:
            grid = (self.mat == k)
        grid = grid.astype(int)

        all_coordinates = []
        for i in range(self.size):
            all_coordinates.append(list(zip([i] * 4, range(4))))
        all_coordinates = [item for sublist in all_coordinates for item in sublist]

        paths = {}
        for i in all_coordinates:
            for j in all_coordinates:
                paths[(i, j)] = minDistance(grid, i[0], i[1], j[0], j[1])

        b = [key for key in list(paths.keys()) if paths[key] == 1000]

        for key in b:
            if paths[(key[1], key[0])] != 1000:
                paths[(key[0], key[1])] = paths[(key[1], key[0])]
        return paths

    def shortest_path(self,sourcex, sourcey, destx, desty, var='biod'):
        '''
        Returns the shortest path between two cells in adjacency matrix associated to graph.
        Default graph considered is biodiversity habitat, can be fuel, or empty.

        :param sourcex: int
        :param sourcey: int
        :param destx: int
        :param desty: int
        :param var: str - "biod","fuel","zeros"
        :return: float
        '''
        if var == "biod":
            grid = (self.mat >= 1)
        elif var == 'fuel':
            grid = (self.mat == 2)
        elif var == 'zero':
            grid = (self.mat == 0)
        grid = grid.astype(int)
        a = minDistance(grid, sourcex, sourcey, destx, desty)
        b = minDistance(grid, destx, desty, sourcey, sourcex)
        if a != b:
            return min(a, b)
        else:
            return a

    def IIC(self, var="biod", nc=True):
        '''
        Integral Index of Connectivity (Pascual Horta, Suara, 2006) with non communicating cells having a finite value
        for shortest path.

        :param var: str - "biod" or "fuel"
        :param nc: bool
        :return: float
        '''
        if nc:
            if var == "biod":
                paths = self.test_algo(1)
            elif var == "fuel":
                paths = self.test_algo(2)
            invert = [1 / (1 + x) for x in list(paths.values())]
            # replace 1/1001 or 1/101 by 0
            iic = sum(invert) / (self.size ** 4)
        else:
            if var == "biod":
                paths = self.test_algo(1)
            elif var == "fuel":
                paths = self.test_algo(2)
            screen1 = [x for x in list(paths.values()) if x < 1000]
            invert = [1 / (1 + x) for x in screen1]
            # replace 1/1001 or 1/101 by 0
            iic = sum(invert) / (self.size ** 4)
        return iic

    def CPL(self, var="biod", nc=True):
        '''
        Characteristic path length associated to graph

        :param var: str - "biod" or "fuel"
        :param nc: bool
        :return: float
        '''
        if nc:
            if var == "biod":
                paths = self.test_algo(1)
            elif var == "fuel":
                paths = self.test_algo(2)
        # need to remove all pairs such that (i=j) from paths
            paths2 = [paths[k] for k in list(paths.keys()) if k[0] != k[1]]
            return stats.mean(paths2)
        else:
            if var == "biod":
                paths = self.test_algo(1)
            elif var == "fuel":
                paths = self.test_algo(2)
            # need to remove all pairs such that (i=j) from paths
            paths2 = [paths[k] for k in list(paths.keys()) if k[0] != k[1]]
            paths2 = [x for x in paths2 if x < 1000]
            try:
                result = stats.mean(paths2)
            except stats.StatisticsError:
                result = float("nan")
            return result

    def diameter(self, var="biod"):
        '''
        Diameter (longest path) associated to graph

        :param var: str - "biod" or "fuel"
        :return: float
        '''
        if var == "zero":
            paths = self.test_algo(0)
        if var == "biod":
            paths = self.test_algo(1)
        elif var == "fuel":
            paths = self.test_algo(2)
        paths2 = [paths[x] for x in paths.keys() if paths[x] < 1000]
        return max(paths2)

    def landscape_shape_index(self, var="biod"):
        '''
        Landscape shape index for raster data, see Fragstat documentation (p. 115)

        :param var: str - "biod" or "fuel"
        :return: float
        '''
        if var == 'biod':
            if self.perimeter("biod" )== 0:
                return 0
            else :
                return 0.25 * self.perimeter(var="biod")/sum(self.shape >= 1)
        elif var == 'fuel':
            if self.perimeter("fuel") == 0:
                return 0
            else:
                return 0.25 * self.perimeter(var="fuel")/sum(self.shape == 2)


    def diagonal(self, direction='NO-SE', var="biod"):
        '''
        Shortest diagonal path

        :param direction: str - "SO-NE"
        :param var: str - "biod" or "fuel"
        :return:
        '''
        if var == "biod":
            paths = self.test_algo(1)
        elif var == "fuel":
            paths = self.test_algo(2)

        if direction == 'NO-SE':
            return paths[((0, 0), (self.size-1, self.size-1))]
        elif direction == 'SO-NE':
            return paths[((self.size-1, 0), (0, self.size-1))]
    # Description
    def number_empty_neighbors(self, i, j, var="biod"):
        '''
        Number of empty neighbor in 4 directions

        :param i: int
            row
        :param j: int
            col
        :param var: str - "biod" or "fuel"
        :return: int
        '''
        size = self.size - 1
        if var == "biod":
            threshold = 1
        elif var == "fuel":
            threshold = 2

        if i > 0:
            up = (self.mat[i-1][j] < threshold)
        elif i == 0:
            up = 1

        if i < size:
            down = (self.mat[i+1][j] < threshold)
        elif i == size:
            down = 1

        if j > 0:
            left = (self.mat[i][j-1] < threshold)
        elif j == 0:
            left = 1

        if j < size:
            right = (self.mat[i][j+1] < threshold)
        elif j == size:
            right = 1

        return sum([up, down, right, left])

    def number_zero_adjacent(self, i, j, var="biod"):
        '''
        Number of empty neighbors in 8 directions

        :param i: int
            row
        :param j: int
            col
        :param var: str - "biod" or "fuel"
        :return: int
        '''

        size = self.size - 1

        if var == "biod":
            threshold = 1
        elif var == "fuel":
            threshold = 2

        if i < 0 or j < 0 or i > size or j > size:
            return 0

        if i > 0:
            up = int(self.mat[i-1][j] < threshold)
        elif i <= 0:
            up = 1

        if i < size:
            down = int(self.mat[i+1][j] < threshold)
        elif i >= size:
            down = 1
        if j > 0:
            left = int(self.mat[i][j-1] < threshold)
        elif j <= 0:
            left = 1
        if j < size:
            right = int(self.mat[i][j+1] < threshold)
        elif j >= size:
            right = 1

        if i > 0 and j < size:
            up_right = int(self.mat[i-1][j+1] < threshold)
        else:
            up_right = 1

        if i > 0 and j > 0:
            up_left = int(self.mat[i-1][j-1] < threshold)
        else:
            up_left = 1

        if i < size and j > 0:
            down_left = int(self.mat[i+1][j-1] < threshold)
        else:
            down_left = 1

        if i < size and j < size:
            down_right = int(self.mat[i+1][j+1] < threshold)
        else:
            down_right = 1

        sum_neighbor = sum([up, down, right, left, up_right, up_left, down_right, down_left])
        return sum_neighbor*(self.mat[i][j] >= threshold)

    def components(self, var="biod"):
        '''
        Returns the number of components/subgraphs in graph of specific variable

        :param var: str - "biod" or "fuel"
        :return: int
        '''
        if var == "biod":
            check = (self.mat >= 1).astype(int)
        elif var == "fuel":
            check = (self.mat >= 2).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)

        return g.countIslands()[0]

    def components_area(self, var='biod'):
        '''
        Returns the area of components/subgraphs associated to specific variable

        :param var: str - "biod" or "fuel"
        :return: int or float
        '''
        if var == "biod":
            check = (self.mat >= 1).astype(int)
        elif var=="fuel":
            check = (self.mat >= 2).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)
        result = g.countIslands()[1]
        if len(result) == 1:
            result = int(result[0])
        elif len(result) == 0:
            result = float("nan")
        return result

    def components_area_max(self, var='biod'):
        '''
        Returns area of largest component/ number of nodes of largest subgraph associated with variable

        :param var: str - "biod" or "fuel"
        :return: int
        '''
        if var == "biod":
            check = (self.mat >= 1).astype(int)
        elif var == "fuel":
            check = (self.mat >=2).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)

        result = g.countIslands()[1]
        if len(result) == 1:
            result = int(result[0])
        elif len(result) == 0:
            result = float("nan")
        else:
            result = max(result)
        return result

    def components_perimeter(self, var="biod", components_var=False):
        '''
        Returns perimeter of components in landscape associated with variable

        :param var: str - "biod" or "fuel"
        :param components_var: bool - optional, enable returning the number of components
        :return: int
        '''
        vars = var
        if var == "biod":
            check = (self.mat >= 1).astype(int)
        elif var == "fuel":
            check = (self.mat >= 2).astype(int)

        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)
        first_nodes = g.countIslands()[2]

        if len(first_nodes) <= 1:
            if var == "biod":
                return self.perimeter(var="biod")
            elif var == 'fuel':
                return self.perimeter(var="fuel")
        else:
            components = []
            perimeter_components = []
            for nodes in first_nodes:
                queue = []
                queue.append(nodes)
                visited = np.ones((4, 4), dtype=bool) == False
                patches = []
                cols = [0, 0, -1, 1, -1, -1, 1, 1]
                rows = [-1, 1, 0, 0, 1, -1, 1, -1]
                while len(queue) > 0:
                    node = queue.pop(0)

                    for i, j in list(zip(rows, cols)):
                        if isValid(node[0] + i, node[1] + j, check, visited):
                            queue.append((node[0] + i, node[1] + j))
                            visited[node[0] + i, node[1] + j] = True
                            patches.append((node[0] + i, node[1] + j))
                components.append(patches)
                perimeter = 0
                for patch in patches:
                    perimeter += check[patch[0], patch[1]] * self.number_empty_neighbors(patch[0], patch[1], var=vars)

                perimeter_components.append(perimeter)
            if components_var:
                return perimeter_components, components
            else:
                return perimeter_components

    def components_shape_index(self, var="biod"):
        '''
        Analogous to Landscape Shape Index at the largest component scale

        :param var: str - "biod" or "fuel"
        :return: float
        '''
        if var == "biod":
            if self.components() > 1:
                area = max(self.components_area(var="biod"))
                candidate = self.components_area(var="biod").index(area)
                perimeter = self.components_perimeter(var="biod")[candidate]
            else:
                area = self.components_area(var="biod")
                perimeter = self.components_perimeter(var="biod")

        elif var == "fuel":
            if self.components("fuel") > 1:
                area = max(self.components_area(var="fuel"))
                candidate = self.components_area(var="fuel").index(area)
                perimeter = self.components_perimeter(var="fuel")[candidate]
            else:
                area = self.components_area(var="fuel")
                perimeter = self.components_perimeter(var="fuel")
        return 0.25*perimeter/math.sqrt(area)

        # Overall graph

    def node_degree(self, i, j, var="biod"):
        '''
        Returns degree of node (i,j)

        :param i: int
            row
        :param j: int
            column
        :param var: str - "biod" or "fuel"
        :return: int
        '''
        if i < 0 or j < 0 or j > self.size - 1 or i > self.size - 1:
            return 0
        if var == "biod":
            if self.mat[i, j] == 0:
                return 0
            else:
                return 8 - self.number_zero_adjacent(i, j, var='biod')
        elif var == "fuel":
            if self.mat[i, j] < 2:
                return 0
            else:
                return 8 - self.number_zero_adjacent(i, j, var='fuel')

    def connectivity_correlation(self, var="biod", option_plot="No"):
        '''
        Returns Pearson coefficient of correlation between node degree and neighbors' degrees.

        :param var: str - "biod" or "fuel"
        :param option_plot: str - optional, returns a graph associated to measure
        :return: float
        '''
        store_node = []
        store_neighbors = []
        for i in range(self.size):
            for j in range(self.size):
                first_coord = [-1, 1, 0, 0, -1, -1, 1, 1]
                second_coord = [0, 0, -1, 1, -1, 1, -1, 1]
                store = []
                if var == "biod":
                    for nb in range(8):
                        store.append(self.node_degree(i + first_coord[nb], j + second_coord[nb], var="biod"))
                    store_node.append(self.node_degree(i, j, var="biod"))
                elif var == "fuel":
                    for nb in range(8):
                        store.append(self.node_degree(i + first_coord[nb], j + second_coord[nb], var="fuel"))
                    store_node.append(self.node_degree(i, j, var="fuel"))

                positive_x = [x for x in store if x != 0]
                try:
                    store_neighbors.append(sum(positive_x)/len(positive_x))
                except stats.StatisticsError:
                    store_neighbors.append(0)
        if len(set(np.array(store_node))) == 1 or len(set(np.array(store_neighbors))) == 1:
            coefficient = np.array([[float('nan'), float('nan')],
                                   [float('nan'), float('nan')]])
        else:
            coefficient = np.corrcoef(np.array(store_node), np.array(store_neighbors))
        if option_plot == "Yes":
            plt.scatter(store_node, store_neighbors)
            plt.ylabel("Average # of edges for 8 neighbors ")
            plt.xlabel("Number of edges of node")
            plt.show()
        else:
            return coefficient[0,1]

    def perimeter(self, var="biod"):
        '''
        Returns perimeter of overall landscape depending on variable

        :param var: str - "biod" or "fuel"
        :return: int
        '''
        perimeter = 0
        if var == "biod":
            for j in range(len(self.mat[0])):
                for i in range(len(self.mat[0])):
                    perimeter += (self.mat[i][j] > 0)*self.number_empty_neighbors(i, j)
        elif var == "fuel":
            for j in range(len(self.mat[0])):
                for i in range(len(self.mat[0])):
                    perimeter += (self.mat[i][j] > 1) * self.number_empty_neighbors(i, j, var='fuel')
        return perimeter

    def perimeter_to_area(self, var='biod'):
        '''
        Returns ratio of perimeter to area on the whole landscape

        :param var: str - "biod" or "fuel"

        :return: float
        '''
        if var == 'biod':
            return self.perimeter(var="biod")/math.sqrt(sum(sum(self.mat >= 1)))
        elif var == 'fuel':
            return self.perimeter(var="fuel")/math.sqrt(sum(sum(self.mat == 2)))

######################################################
def biod(shape):
    '''
    Returns biodiversity score defined in working paper
    :param shape: np.array
        land
    :return: int
        score
    '''
    y = (shape >= 1).dot(utilities.connectivity_mat()).dot(np.transpose((shape >= 1)))
    return y
#biod(shape)

def biod_relative(shape):
    """
    Normalized biodiversity score as a share of largest biodiversity score
    :param shape: np.array
        land
    :return: int
        score
    """
    y = biod(shape)/biod(np.repeat(1,params.size**2))
    return y
#biod_relative(shape)

def fuel(shape):
    """
    Returns fuel connectivity score defined in working paper
    :param shape: np.array
        land
    :return: int
        score
    """
    y = (shape >= 2).dot(utilities.connectivity_mat()).dot(np.transpose((shape >= 2)))
    return y
#fuel(shape)

def fuel_relative(shape):
    """
    Normalized fuel score as a share of maximum score possible
    :param shape: np.array
        land
    :return: int
        score
    """
    y = fuel(shape)/fuel(np.repeat(2,params.size**2))
    return y
#fuel_relative(shape)

def node_biod(shape):
    '''
    Returns number of nodes in the biodiversity graph
    :param shape: np.array
        land
    :return: int
        nb nodes
    '''
    y = sum(shape>0)
    return y
#node_biod(shape)

def node_fuel(shape):
    '''
    Returns number of nodes in the fuel graph
    :param shape: np.array
        land
    :return: int
        nb nodes
    '''
    y = sum(shape==2)
    return y
#node_fuel(shape)

def habitat_ratio(shape):
    '''
    Indicator of biodiversity relative to overall landscape area
    :param shape: np.array
        land
    :return: int
        score
    '''
    y = sum(shape>0)/len(shape)
    return y
#habitat_ratio(shape)

def fuel_ratio(shape):
    """
    Indicator of fuel relative to overall landscape area
    :param shape: np.array
        land
    :return: int
        score
    """
    y = node_fuel(shape)/len(shape)
    return y
#fuel_ratio(shape)

def test_algo(shape, k):
    """
    Fundamental algorithm computing the minimum distance between nodes of graphs (biodiversity or fuel)
    :param shape: np.array
        shape
    :param k: int
        threshold for fuel or biodiversity habitat
    :return: dict
        nodes coordinates and distances
    """
    if k >= 0:
        grid = utilities.to_matrix(shape>=k)
    elif k==0:
        grid = utilities.to_matrix(shape == k)
    grid = grid.astype(int)

    all_coordinates = []
    for i in range(int(math.sqrt(len(shape)))):
        all_coordinates.append(list(zip([i] * params.size, range(params.size))))
    all_coordinates = [item for sublist in all_coordinates for item in sublist]

    paths = {}
    for i in all_coordinates:
        for j in all_coordinates:
            paths[(i, j)] = minDistance(grid, i[0], i[1], j[0], j[1])

    b = [key for key in list(paths.keys()) if paths[key] == 1000]

    for key in b:
        if paths[(key[1], key[0])] != 1000:
            paths[(key[0], key[1])] = paths[(key[1], key[0])]
    return paths
#test_algo(shape,0)
#test_algo(shape,1)
#test_algo(shape,2)

def number_empty_neighbors(shape, i, j, var="biod"):
    '''
    Number of empty neighbor in 4 directions

    :param i: int
        row
    :param j: int
        col
    :param var: str - "biod" or "fuel"
    :return: int
    '''
    size = params.size - 1
    if var == "biod":
        threshold = 1
    elif var == "fuel":
        threshold = 2

    if i > 0:
        up = (utilities.to_matrix(shape)[i - 1][j] < threshold)
    elif i == 0:
        up = 1

    if i < size:
        down = (utilities.to_matrix(shape)[i + 1][j] < threshold)
    elif i == size:
        down = 1

    if j > 0:
        left = (utilities.to_matrix(shape)[i][j - 1] < threshold)
    elif j == 0:
        left = 1

    if j < size:
        right = (utilities.to_matrix(shape)[i][j + 1] < threshold)
    elif j == size:
        right = 1

    return sum([up, down, right, left])
#number_empty_neighbors(shape, 1, 2)

def perimeter(shape, var="biod"):
    """
    Perimeter of a landscape's biodiversity habitat or fuel
    :param shape: np.array
        land
    :param var: str
        "biod" or "fuel"
    :return: int
        Perimeter
    """
    perimeter = 0
    if var == "biod":
        for j in range(len(utilities.to_matrix(shape)[0])):
            for i in range(len(utilities.to_matrix(shape)[0])):
                perimeter += (utilities.to_matrix(shape)[i][j] > 0)*number_empty_neighbors(shape, i, j)
    elif var == "fuel":
        for j in range(len(utilities.to_matrix(shape)[0])):
            for i in range(len(utilities.to_matrix(shape)[0])):
                perimeter += (utilities.to_matrix(shape)[i][j] > 1) * number_empty_neighbors(shape, i, j, var='fuel')
    return perimeter
#perimeter(shape)

def shortest_path(shape,sourcex, sourcey, destx, desty, var='biod'):
    '''
    Returns the shortest path between two cells in adjacency matrix associated to graph.
    Default graph considered is biodiversity habitat, can be fuel, or empty.

    :param sourcex: int
    :param sourcey: int
    :param destx: int
    :param desty: int
    :param var: str
    "biod","fuel","zeros"
    :return: float
    '''
    if var == "biod":
        grid = (utilities.to_matrix(shape) >= 1)
    elif var == 'fuel':
        grid = (utilities.to_matrix(shape) == 2)
    elif var == 'zero':
        grid = (utilities.to_matrix(shape) == 0)
    grid = grid.astype(int)
    a = minDistance(grid, sourcex, sourcey, destx, desty)
    b = minDistance(grid, destx, desty, sourcey, sourcex)
    if a != b:
        return min(a, b)
    else:
        return a
#shortest_path(shape, 0,0,3,2)

def IIC(shape, var="biod", nc=True):
    '''
    Integral Index of Connectivity (Pascual Horta, Suara, 2006) with non communicating cells having a finite value
    for shortest path.

    :param var: str
        Alternative : 'fuel'
    :param nc: bool
    :return: float
    '''
    if nc:
        if var == "biod":
            paths = test_algo(shape,1)
        elif var == "fuel":
            paths = test_algo(shape,2)
        invert = [1 / (1 + x) for x in list(paths.values())]
            # replace 1/1001 or 1/101 by 0
        iic = sum(invert) / (params.size ** 4)
    else:
        if var == "biod":
            paths = test_algo(shape,1)
        elif var == "fuel":
            paths = test_algo(shape,2)
        screen1 = [x for x in list(paths.values()) if x < 1000]
        invert = [1 / (1 + x) for x in screen1]
            # replace 1/1001 or 1/101 by 0
        iic = sum(invert) / (params.size ** 4)
    return iic
#IIC(shape)

def PC(shape):
    """
    Implements Probability of Connectivity following Saura and Pascual Horta, 2007
    :param shape:
    :return:
    """
    paths = test_algo(shape, 1)
    dist = [x for x in list(paths.values())]
    prob = [math.exp(-x) for x in dist]
    return sum(prob)/(params.size**4)


def CPL(shape, var="biod", nc=True):
    '''
    Characteristic path length associated to graph

    :param var: str- "biod" or "fuel"
    :param nc: bool
    :return: float
    '''
    if nc:
        if var == "biod":
            paths = test_algo(shape, 1)
        elif var == "fuel":
            paths = test_algo(shape, 2)
            # need to remove all pairs such that (i=j) from paths
        paths2 = [paths[k] for k in list(paths.keys()) if k[0] != k[1]]
        return stats.mean(paths2)
    else:
        if var == "biod":
            paths = test_algo(shape, 1)
        elif var == "fuel":
            paths = test_algo(shape, 2)
            # need to remove all pairs such that (i=j) from paths
        paths2 = [paths[k] for k in list(paths.keys()) if k[0] != k[1]]
        paths2 = [x for x in paths2 if x < 1000]
        try:
            result = stats.mean(paths2)
        except stats.StatisticsError:
            result = float("nan")
        return result
#CPL(shape)

def diameter(shape, var="biod"):
    '''
    Diameter (longest path) associated to graph

    :param var: str - "biod" or "fuel"
    :return: float
    '''
    if var == "zero":
        paths = test_algo(shape, 0)
    if var == "biod":
        paths = test_algo(shape, 1)
    elif var == "fuel":
        paths = test_algo(shape, 2)
    paths2 = [paths[x] for x in paths.keys() if paths[x] < 1000]
    return max(paths2)
#diameter(shape)

def landscape_shape_index(shape, var="biod"):
    '''
    Landscape shape index for raster data, see Fragstat documentation (p. 115)

    :param var: str - "biod" or "fuel"
    :return: float
    '''
    if var == 'biod':
        if perimeter(shape, "biod" )== 0:
            return 0
        else :
            return 0.25 * perimeter(shape, var="biod")/math.sqrt(sum(shape >= 1))
    elif var == 'fuel':
        if perimeter(shape, "fuel") == 0:
            return 0
        else:
            return 0.25 * perimeter(shape, var="fuel")/math.sqrt(sum(shape == 2))

    # Description
    # Utility functions : empty neighbors and zero adjacent
#landscape_shape_index(shape)

def number_zero_adjacent(shape, i, j, var="biod"):
    '''
    Number of empty neighbors in 8 directions

    :param i: int
        row
    :param j: int
        col
    :param var: str - "biod" or "fuel"
    :return: int
    '''

    size = params.size - 1

    if var == "biod":
        threshold = 1
    elif var == "fuel":
        threshold = 2

    if i < 0 or j < 0 or i > size or j > size:
        return 0

    if i > 0:
        up = int(utilities.to_matrix(shape)[i - 1][j] < threshold)
    elif i <= 0:
        up = 1

    if i < size:
        down = int(utilities.to_matrix(shape)[i + 1][j] < threshold)
    elif i >= size:
        down = 1
    if j > 0:
        left = int(utilities.to_matrix(shape)[i][j - 1] < threshold)
    elif j <= 0:
        left = 1
    if j < size:
        right = int(utilities.to_matrix(shape)[i][j + 1] < threshold)
    elif j >= size:
        right = 1

    if i > 0 and j < size:
        up_right = int(utilities.to_matrix(shape)[i - 1][j + 1] < threshold)
    else:
        up_right = 1

    if i > 0 and j > 0:
        up_left = int(utilities.to_matrix(shape)[i - 1][j - 1] < threshold)
    else:
        up_left = 1

    if i < size and j > 0:
        down_left = int(utilities.to_matrix(shape)[i + 1][j - 1] < threshold)
    else:
        down_left = 1

    if i < size and j < size:
        down_right = int(utilities.to_matrix(shape)[i + 1][j + 1] < threshold)
    else:
        down_right = 1

    sum_neighbor = sum([up, down, right, left, up_right, up_left, down_right, down_left])
    return sum_neighbor * (utilities.to_matrix(shape)[i][j] >= threshold)
#number_zero_adjacent(shape, 1, 3)

def components(shape, var="biod"):
    '''
    Returns number, area and initial coordinates of components/subgraphs, depending on the variable

    :param var: str - "biod" or "fuel"
    :return: list
    '''
    if var == "biod":
        check = (utilities.to_matrix(shape) >= 1).astype(int)
    elif var == "fuel":
        check = (utilities.to_matrix(shape) >= 2).astype(int)
    graph = list([list(x) for x in check])
    row = len(graph)
    col = len(graph[0])

    g = Graph(row, col, graph)
    return g.countIslands()
#components(shape)

def components_perimeter(shape, var="biod", components_var=False):
    '''
    Returns the perimeters associated with components/subgraphs depending on the variable considered

    :param var: str - "biod" or "fuel"
    :param components_var: bool - optional, returns components as well
    :return: list int
    '''
    vars = var
    first_nodes = 0
    if var == "biod":
        check = (utilities.to_matrix(shape) >= 1).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)

        first_nodes = g.countIslands()[2]
        check = (utilities.to_matrix(shape) >= 1).astype(int)
    elif var == "fuel":
        check = (utilities.to_matrix(shape) >= 2).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)

        first_nodes = g.countIslands()[2]
        check = (utilities.to_matrix(shape) == 2).astype(int)
    if len(first_nodes) <= 1:
        perimeter = 0
        if var == "biod":
            for j in range(len(utilities.to_matrix(shape)[0])):
                for i in range(len(utilities.to_matrix(shape)[0])):
                    perimeter += (utilities.to_matrix(shape)[i][j] > 0) * number_empty_neighbors(shape, i, j)
            return perimeter
        elif var == 'fuel':
            for j in range(len(utilities.to_matrix(shape)[0])):
                for i in range(len(utilities.to_matrix(shape)[0])):
                    perimeter += (utilities.to_matrix(shape)[i][j] > 1) * number_empty_neighbors(shape, i, j,
                                                                                                 var='fuel')
            return perimeter
    else:
        components = []
        perimeter_components = []
        for nodes in first_nodes:
            queue = []
            queue.append(nodes)
            visited = np.zeros((4, 4), dtype=bool)
            patches = []
            cols = [0, 0, -1, 1, -1, -1, 1, 1]
            rows = [-1, 1, 0, 0, 1, -1, 1, -1]
            while len(queue) > 0:
                node = queue.pop(0)

                for i, j in list(zip(rows, cols)):
                    if isValid(node[0] + i, node[1] + j, check, visited):
                        queue.append((node[0] + i, node[1] + j))
                        visited[node[0] + i, node[1] + j] = True
                        patches.append((node[0] + i, node[1] + j))
            components.append(patches)
            perimeter = 0
            for patch in patches:
                perimeter += check[patch[0], patch[1]] * number_empty_neighbors(shape, patch[0], patch[1], var=vars)

            perimeter_components.append(perimeter)
        if components_var:
            return perimeter_components, components
        else:
            return perimeter_components
#components_perimeter(shape)

def components_shape_index(shape, var="biod"):
    '''
    Analogous to Landscape Shape Index at the largest component scale

    :param var: str - "biod" or "fuel"
    :return: float
    '''

    if var == "biod":
        y = components(shape, var = 'biod')
        if y[0] > 1:
            area = max(y[1])
            candidate = y[1].index(area)
            try:
                perimeter = components_perimeter(shape, var="biod")[candidate]
            except IndexError:
                perimeter = components_perimeter(shape, var="biod")
        elif y[0]==0:
            return 0
        else:
            area = y[1][0]
            perimeter = components_perimeter(shape, var="biod")

    elif var == "fuel":
        y = components(shape, var='fuel')
        if y[0] > 1:
            area = max(y[1])
            candidate = y[1].index(area)
            try:
                perimeter = components_perimeter(shape, var="fuel")[candidate]
            except IndexError:
                perimeter = components_perimeter(shape, var="fuel")

        elif y[0]==0:
            return 0
        else:
            area = y[1][0]
            perimeter = components_perimeter(shape, var="fuel")
    return 0.25 * perimeter/math.sqrt(area)

        # Overall graph
#components_shape_index(shape)

def node_degree(shape, i, j, var="biod"):
    '''
    Returns degree of node (i,j)

    :param i: int
        row
    :param j: int
        column
    :param var: str - "biod" or "fuel"
    :return: int
    '''
    if i < 0 or j < 0 or j > params.size-1 or i > params.size-1:
        return 0
    if var == "biod":
        if utilities.to_matrix(shape)[i, j] == 0:
            return 0
        else:
            return 8 - number_zero_adjacent(shape, i, j, var='biod')
    elif var == "fuel":
        if utilities.to_matrix(shape)[i, j] < 2:
            return 0
        else:
            return 8 - number_zero_adjacent(shape, i, j, var='fuel')
#node_degree(shape, 1, 2)


def equivalent(shape):
    m2 = np.rot90(utilities.to_matrix(shape))
    m3 = np.rot90(m2)
    m4 = np.rot90(m3)
    if params.size == 5:
        # Horizontal symmetry
        hor_sym = [utilities.to_matrix(shape)[4], utilities.to_matrix(shape)[3], utilities.to_matrix(shape)[2], utilities.to_matrix(shape)[1], utilities.to_matrix(shape)[0]]
    elif params.size == 4:
        hor_sym = [utilities.to_matrix(shape)[3], utilities.to_matrix(shape)[2], utilities.to_matrix(shape)[1], utilities.to_matrix(shape)[0]]
    elif params.size == 3:
        hor_sym = [utilities.to_matrix(shape)[2], utilities.to_matrix(shape)[1], utilities.to_matrix(shape)[0]]
        # 3 clockwise transformations
    else:
        ValueError("Size is not supported")
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
#equivalent(shape)

def convergence(successions):
    """
    Returns the time of convergence and length of cycle if applicable

    :param successions: data series
    :return: a list of tuples (t_convergence, cycle_period)
    """
    lands = [Land(shape=np.array(ast.literal_eval(x))) for x in successions]
    horizon = list(range(len(successions)))
    horizon.remove(0)
    for t in horizon:
        for i in range(1,t):
            if len(lands[t].equivalent().intersection(lands[t - i].equivalent())) > 0:
                break
        else:
            continue
        break
    convergence_time = t
    initial = t-i
    cycle = i

    for period in range(initial,horizon[-1]-cycle,cycle):
        if len(lands[period].equivalent().intersection(lands[period+cycle].equivalent()))==0:
            print('Not a convergence')
            #return float('nan'), float('nan')
            break

    else:
        #print('Convergence in ' + str(cycle))
        return convergence_time, cycle
#successions = [str(shape)]*10

def data_fix(data_path):
    """
    Small function to fix some stats that were wrong before
    :param data_path: str
    :return: .csv
        statistics
    """
    data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/sample_data/"+ data_path)

    values = list(data.value)
    values_unique = set(values)

    land_CPL_biod = {}
    land_CPL_fuel = {}
    land_shape_index_biod = {}
    land_shape_index_fuel = {}

    for shape in values_unique:
        land = np.array(ast.literal_eval(shape))
        land_CPL_biod[shape] = CPL(land, nc=False)
        land_CPL_fuel[shape] = CPL(land, var="fuel", nc=False)
        land_shape_index_biod[shape] = landscape_shape_index(land)
        land_shape_index_fuel[shape] = landscape_shape_index(land, var="fuel")

    data["land_CPL_biod"] = [land_CPL_biod.get(n, n) for n in values]
    data["land_CPL_fuel"] = [land_CPL_fuel.get(n, n) for n in values]
    data["land_shape_index_biod"] = [land_shape_index_biod.get(n,n) for n in values]
    data['land_shape_index_fuel'] = [land_shape_index_fuel.get(n,n) for n in values]

    col = list(data.columns).index("value")
    for x in range(20, len(data) + 1, 20):
        succ = data.iloc[x - 20:x, col]
        conv = convergence(succ)
        data.iloc[x - 20:x, list(data.columns).index('convergence_time')] = conv[0]
        data.iloc[x - 20:x, list(data.columns).index('convergence_period')] = conv[1]

    data.to_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/sample_data/" + data_path)

def statistics_dataframe(data_path, output_=False):
    """
    Computes a wide array of statistics for a whole dataset
    :param data_path: str
        Path to data
    :param output_:
    :return: .csv
    """
    index = re.findall(r'\d+', data_path)
    #budget_here = index[2]
    if params.size==3:
        path = params.path+ "successions_matched/"
    elif params.size==4:
        path = params.path + "successions_matched/"

    # Load data
    data = pd.read_csv(path+data_path).drop(columns=['Unnamed: 0'])

    # Keep 20 steps
    data = data.iloc[:, 0:20]
    names = list(range(20))
    data.columns = names

    # keep track of landscapes with index
    data["index"] = list(range(len(data)))
    checker = data.melt(id_vars=["index"])
    checker["variable"] = checker['variable'].astype(int)
    checker = checker.sort_values(["index", "variable"])
    checker = checker.reset_index(drop=True)

    # List values of landscapes (np.array) and only keep unique, to minimize the number of operations
    values = list(checker.value)
    values_unique = list(set(values))

    #region Setting up of data storage dictionnaries
    nodes_biod = {}
    nodes_fuel = {}

    biod_score = {}
    fuel_score = {}

    land_perimeter_biod = {}
    land_perimeter_fuel = {}

    land_habitat_area_ratio = {}
    land_fuel_area_ratio = {}

    land_shape_index_biod = {}
    land_shape_index_fuel = {}

    land_diameter_biod = {}
    land_diameter_fuel = {}

    land_IIC_nc_biod = {}
    land_IIC_nc_fuel = {}

    land_IIC_biod = {}
    land_IIC_fuel = {}

    land_CPL_nc_biod = {}
    land_CPL_nc_fuel = {}

    land_CPL_biod = {}
    land_CPL_fuel = {}

    components_number_biod = {}
    components_number_fuel = {}

    components_area_biod_max = {}
    components_area_fuel_max = {}

    components_shape_index_biod = {}
    components_shape_index_fuel = {}

    land_PC = {}
    # endregion

    #region Compute statistics
    for shape in values_unique:
        land = np.array(ast.literal_eval(shape))

        nodes_biod[shape] = node_biod(land)
        nodes_fuel[shape] = node_fuel(land)

        biod_score[shape] = biod(land)
        fuel_score[shape] = fuel(land)

        land_perimeter_biod[shape] = perimeter(land)
        land_perimeter_fuel[shape] = perimeter(land, var="fuel")

        land_habitat_area_ratio[shape] = habitat_ratio(land)
        land_fuel_area_ratio[shape] = fuel_ratio(land)

        land_shape_index_biod[shape] = landscape_shape_index(land)
        land_shape_index_fuel[shape] = landscape_shape_index(land, var="fuel")

        land_diameter_biod[shape] = diameter(land)
        land_diameter_fuel[shape] = diameter(land, var='fuel')

        land_IIC_nc_biod[shape] = IIC(land)
        land_IIC_nc_fuel[shape] = IIC(land, var="fuel")

        land_IIC_biod[shape] = IIC(land, nc=False)
        land_IIC_fuel[shape] = IIC(land, var="fuel", nc=False)

        land_CPL_nc_biod[shape] = CPL(land)
        land_CPL_nc_fuel[shape] = CPL(land, var='fuel')

        land_CPL_biod[shape] = CPL(land, nc=False)
        land_CPL_fuel[shape] = CPL(land, var="fuel", nc=False)

        components_number_biod[shape] = components(land)[0]
        components_number_fuel[shape] = components(land, var="fuel")[0]

        try:
            components_area_biod_max[shape] = max(components(land)[1])
        except ValueError:
            components_area_biod_max[shape] = 0
        try:
            components_area_fuel_max[shape] = max(components(land, var="fuel")[1])
        except ValueError:
            components_area_fuel_max[shape] = 0

        components_shape_index_biod[shape] = components_shape_index(land)
        components_shape_index_fuel[shape] = components_shape_index(land, var="fuel")

        land_PC[shape] = PC(land)
    #endregion

    #region Assign values to dataset
    checker["nodes_biod"] = [nodes_biod.get(n, n) for n in values]
    checker["nodes_fuel"] = [nodes_fuel.get(n, n) for n in values]
    checker["nodes_zero"] = [16 - checker.nodes_biod.loc[x] for x in range(len(checker))]
    checker["score_biod"] = [biod_score.get(n, n) for n in values]
    checker["score_fuel"] = [fuel_score.get(n, n) for n in values]

    checker["land_perimeter_biod"] = [land_perimeter_biod.get(n, n) for n in values]
    checker['land_perimeter_fuel'] = [land_perimeter_fuel.get(n, n) for n in values]
    checker["land_shape_index_biod"] = [land_shape_index_biod.get(n, n) for n in values]
    checker["land_shape_index_fuel"] = [land_shape_index_fuel.get(n, n) for n in values]

    checker["land_diameter_biod"] = [land_diameter_biod.get(n, n) for n in values]
    checker["land_diameter_fuel"] = [land_diameter_fuel.get(n, n) for n in values]

    checker["land_IIC_nc_biod"] = [land_IIC_nc_biod.get(n, n) for n in values]
    checker["land_IIC_nc_fuel"] = [land_IIC_nc_fuel.get(n, n) for n in values]

    checker["land_IIC_biod"] = [land_IIC_biod.get(n, n) for n in values]
    checker["land_IIC_fuel"] = [land_IIC_fuel.get(n, n) for n in values]

    checker["land_CPL_nc_biod"] = [land_CPL_nc_biod.get(n, n) for n in values]
    checker["land_CPL_nc_fuel"] = [land_CPL_nc_fuel.get(n, n) for n in values]

    checker["land_CPL_biod"] = [land_CPL_biod.get(n, n) for n in values]
    checker["land_CPL_fuel"] = [land_CPL_fuel.get(n, n) for n in values]

    checker["components_number_biod"] = [components_number_biod.get(n, n) for n in values]
    checker["components_number_fuel"] = [components_number_fuel.get(n, n) for n in values]

    checker["components_area_max_biod"] = [components_area_biod_max.get(n, n) for n in values]
    checker["components_area_max_fuel"] = [components_area_fuel_max.get(n, n) for n in values]

    checker["components_shape_index_biod"] = [components_shape_index_biod.get(n, n) for n in values]
    checker["components_shape_index_fuel"] = [components_shape_index_fuel.get(n, n) for n in values]

    convergence_time = []
    convergence_period = []
    indices = list(range(0, len(checker) + 1, 20))
    indices.remove(0)
    for x in indices:
        succession = checker.value.loc[x - 20:x - 1]
        result = convergence(succession)
        convergence_time_s = result[0]
        convergence_period_s = result[1]
        convergence_time.extend([convergence_time_s] * 20)
        convergence_period.extend([convergence_period_s] * 20)
    del result

    checker['PC'] = [land_PC.get(n,n) for n in values]

    checker['convergence_time'] = convergence_time
    checker['convergence_period'] = convergence_period

    checker['index_biod'] = [index[-1]] * len(checker)
    checker['budget'] = [params.budget] * len(checker)
    #endregion

    #Save data
    checker.to_csv(params.path_to_save_results +"/budget_" + str(params.budget) + "/statistics/stats_"+data_path)
    if output_ == True:
        return checker.head()

################

def view3D(data,x,y,z,rotation1,rotation2,legend=False):
    '''
    Returns the envelope graph of a function in 3D

    :param data: dataset
    :param x: column
    :param y: column
    :param z: column
    :param rotation1: int
        Vertical angle
    :param rotation2: int
        Horizontal angle
    :param legend: bool
        If True, show color gradient in legend
    :return: plt
    '''
    x_1 = np.linspace(data[x].min(), data[x].max(), len(data[x].unique()))
    y_1 = np.linspace(data[y].min(), data[y].max(), len(data[y].unique()))

    x_2,y_2 = np.meshgrid(x_1,y_1)

    z2 = scipy.interpolate.griddata((data[x],data[y]), data[z],(x_2,y_2), method='cubic')

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    surf = ax.plot_surface(x_2, y_2, z2, rstride=1, cstride=1, cmap=matplotlib.cm.coolwarm,
                       linewidth=0, antialiased=False)
    ax.set_zlim(0,100)

    ax.zaxis.set_major_locator(matplotlib.ticker.LinearLocator(10))
    ax.zaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.02f'))

    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.xlabel(x)
    plt.ylabel(y)

    if legend == True:
        fig.colorbar(surf, shrink=0.5, aspect=10,location='right')
    plt.title('3x3 landscape - value function ')
    ax.view_init(rotation1,rotation2)
    plt.show()


def statistics_dataframe_centrality(data_path):
    """
    Returns statistics for data resulting from centrality algorithm. Different function because slightly different data management approach
    :param data_path: str
        path to data
    :return: .csv
        statistics
    """
    # List values of landscapes (np.array) and only keep unique, to minimize the number of operations
    checker = pd.read_csv(params.path_to_save_results + "/budget_" + str(params.budget) + "/successions_centrality/" + data_path)

    values = list(checker.value)
    values_unique = list(set(values))

    #region Setting up of data storage dictionnaries
    nodes_biod = {}
    nodes_fuel = {}

    biod_score = {}
    fuel_score = {}

    land_perimeter_biod = {}
    land_perimeter_fuel = {}

    land_habitat_area_ratio = {}
    land_fuel_area_ratio = {}

    land_shape_index_biod = {}
    land_shape_index_fuel = {}

    land_diameter_biod = {}
    land_diameter_fuel = {}

    land_IIC_nc_biod = {}
    land_IIC_nc_fuel = {}

    land_IIC_biod = {}
    land_IIC_fuel = {}

    land_CPL_nc_biod = {}
    land_CPL_nc_fuel = {}

    land_CPL_biod = {}
    land_CPL_fuel = {}

    components_number_biod = {}
    components_number_fuel = {}

    components_area_biod_max = {}
    components_area_fuel_max = {}

    components_shape_index_biod = {}
    components_shape_index_fuel = {}

    land_PC = {}
    # endregion

    #region Compute statistics
    for shape in values_unique:
        land = np.array(ast.literal_eval(shape.replace(' ',',')))

        nodes_biod[shape] = node_biod(land)
        nodes_fuel[shape] = node_fuel(land)

        biod_score[shape] = biod(land)
        fuel_score[shape] = fuel(land)

        land_perimeter_biod[shape] = perimeter(land)
        land_perimeter_fuel[shape] = perimeter(land, var="fuel")

        land_habitat_area_ratio[shape] = habitat_ratio(land)
        land_fuel_area_ratio[shape] = fuel_ratio(land)

        land_shape_index_biod[shape] = landscape_shape_index(land)
        land_shape_index_fuel[shape] = landscape_shape_index(land, var="fuel")

        land_diameter_biod[shape] = diameter(land)
        land_diameter_fuel[shape] = diameter(land, var='fuel')

        land_IIC_nc_biod[shape] = IIC(land)
        land_IIC_nc_fuel[shape] = IIC(land, var="fuel")

        land_IIC_biod[shape] = IIC(land, nc=False)
        land_IIC_fuel[shape] = IIC(land, var="fuel", nc=False)

        land_CPL_nc_biod[shape] = CPL(land)
        land_CPL_nc_fuel[shape] = CPL(land, var='fuel')

        land_CPL_biod[shape] = CPL(land, nc=False)
        land_CPL_fuel[shape] = CPL(land, var="fuel", nc=False)

        components_number_biod[shape] = components(land)[0]
        components_number_fuel[shape] = components(land, var="fuel")[0]

        try:
            components_area_biod_max[shape] = max(components(land)[1])
        except ValueError:
            components_area_biod_max[shape] = 0
        try:
            components_area_fuel_max[shape] = max(components(land, var="fuel")[1])
        except ValueError:
            components_area_fuel_max[shape] = 0

        components_shape_index_biod[shape] = components_shape_index(land)
        components_shape_index_fuel[shape] = components_shape_index(land, var="fuel")

        land_PC[shape] = PC(land)
    #endregion

    #region Assign values to dataset
    checker["nodes_biod"] = [nodes_biod.get(n, n) for n in values]
    checker["nodes_fuel"] = [nodes_fuel.get(n, n) for n in values]
    checker["nodes_zero"] = [16 - checker.nodes_biod.loc[x] for x in range(len(checker))]
    checker["score_biod"] = [biod_score.get(n, n) for n in values]
    checker["score_fuel"] = [fuel_score.get(n, n) for n in values]

    checker["land_perimeter_biod"] = [land_perimeter_biod.get(n, n) for n in values]
    checker['land_perimeter_fuel'] = [land_perimeter_fuel.get(n, n) for n in values]
    checker["land_shape_index_biod"] = [land_shape_index_biod.get(n, n) for n in values]
    checker["land_shape_index_fuel"] = [land_shape_index_fuel.get(n, n) for n in values]

    checker["land_diameter_biod"] = [land_diameter_biod.get(n, n) for n in values]
    checker["land_diameter_fuel"] = [land_diameter_fuel.get(n, n) for n in values]

    checker["land_IIC_nc_biod"] = [land_IIC_nc_biod.get(n, n) for n in values]
    checker["land_IIC_nc_fuel"] = [land_IIC_nc_fuel.get(n, n) for n in values]

    checker["land_IIC_biod"] = [land_IIC_biod.get(n, n) for n in values]
    checker["land_IIC_fuel"] = [land_IIC_fuel.get(n, n) for n in values]

    checker["land_CPL_nc_biod"] = [land_CPL_nc_biod.get(n, n) for n in values]
    checker["land_CPL_nc_fuel"] = [land_CPL_nc_fuel.get(n, n) for n in values]

    checker["land_CPL_biod"] = [land_CPL_biod.get(n, n) for n in values]
    checker["land_CPL_fuel"] = [land_CPL_fuel.get(n, n) for n in values]

    checker["components_number_biod"] = [components_number_biod.get(n, n) for n in values]
    checker["components_number_fuel"] = [components_number_fuel.get(n, n) for n in values]

    checker["components_area_max_biod"] = [components_area_biod_max.get(n, n) for n in values]
    checker["components_area_max_fuel"] = [components_area_fuel_max.get(n, n) for n in values]

    checker["components_shape_index_biod"] = [components_shape_index_biod.get(n, n) for n in values]
    checker["components_shape_index_fuel"] = [components_shape_index_fuel.get(n, n) for n in values]

    convergence_time = []
    convergence_period = []
    indices = list(range(0, len(checker) + 1, 20))
    indices.remove(0)
    for x in indices:
        succession = checker.value.loc[x - 20:x - 1]
        succession = [x.replace(' ',',') for x in succession]
        result = convergence(succession)
        convergence_time_s = result[0]
        convergence_period_s = result[1]
        convergence_time.extend([convergence_time_s] * 20)
        convergence_period.extend([convergence_period_s] * 20)
    del result

    checker['PC'] = [land_PC.get(n,n) for n in values]

    checker['convergence_time'] = convergence_time
    checker['convergence_period'] = convergence_period

    #endregion

    #Save data
    checker.to_csv(params.path_to_save_results +"/budget_" + str(params.budget) + "/stat_centrality/stats_centrality_"+data_path)
