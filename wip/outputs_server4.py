import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import ast
import statistics
import random
import os
import re
import csv
import math
import statistics as stats


import utilities_server as utilities
import params_server4 as params


class Graph:

    def __init__(self, row, col, g):
        self.ROW = row
        self.COL = col
        self.graph = g

    # A function to check if a given cell
    # (row, col) can be included in DFS
    def isSafe(self, i, j, visited):
        # row number is in range, column number
        # is in range and value is 1
        # and not yet visited
        return (i >= 0 and i < self.ROW and
                j >= 0 and j < self.COL and
                not visited[i][j] and self.graph[i][j])

    # A utility function to do DFS for a 2D
    # boolean matrix. It only considers
    # the 8 neighbours as adjacent vertices

    def DFS(self, i, j, visited):

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

    # The main function that returns
    # count of islands in a given boolean
    # 2D matrix

    def countIslands(self):
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
                    #area.append(sum([sum(x) for x in visited])-area[-1])
                    area.append(sum([sum(x) for x in visited])-sum(area))
                    zone_l.append((i, j))
        area.pop(0)

        return count, area, zone_l
class QItem:
    def __init__(self, row, col, dist):
        self.row = row
        self.col = col
        self.dist = dist

    def __repr__(self):
        return f"QItem({self.row}, {self.col}, {self.dist})"

class Land:
    def __init__(self, shape):
        self.shape = shape
        self.size = int(math.sqrt(len(shape)))
        self.mat = utilities.to_matrix(shape)
        self.biod = (shape >= 1).dot(utilities.connectivity_mat()).dot(np.transpose((shape >= 1)))
        self.biod_relative = self.biod/(np.repeat(1, len(self.shape)).dot(utilities.connectivity_mat()).dot(np.transpose(np.repeat(1, len(self.shape)))))
        self.fuel = (shape >= 2).dot(utilities.connectivity_mat()).dot(np.transpose((shape >= 2)))
        self.fuel_relative = self.fuel/(np.repeat(1, len(self.shape)).dot(utilities.connectivity_mat()).dot(np.transpose(np.repeat(1, len(self.shape)))))
        self.node_biod = sum(shape > 0)
        self.node_fuel = sum(shape == 2)
        self.zero = sum(shape == 0)
        self.habitat_area_ratio = sum(self.shape >= 1)/(self.size**2)
        self.fuel_area_ratio = sum(self.shape == 2)/(self.size**2)

    def view(self):
        ax = plt.axes()
        ax.set_title('Initial biodiv = '+str(self.biod))
        sns.heatmap(self.mat, ax= ax, linewidth=0.5, cmap="Greens", vmin=0, vmax=2, center=1)
        plt.show()

    def view_succession(self, index):
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
        post = (np.array(self.shape) + 1) * (1 - np.array(presc_burn))
        land_post = np.minimum(post, 2)
        self.fuel_dynamics = land_post

    def succession(self, constraint_biodiv):
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

    def shortest_path(self, sourcex, sourcey, destx, desty):
        grid = (self.mat >= 1)
        grid = grid.astype(int)
        a = minDistance(grid, sourcex, sourcey, destx, desty)
        b = minDistance(grid, destx, desty, sourcey, sourcex)
        if a != b:
            return min(a, b)
        else:
            return a

    def IIC_with_nc(self, var="biod"):
        if var == "biod":
            paths = self.test_algo(1)
        elif var == "fuel":
            paths = self.test_algo(2)
        invert = [1 / (1 + x) for x in list(paths.values())]
        # replace 1/1001 or 1/101 by 0
        iic = sum(invert) / (self.size ** 4)
        return iic

    def IIC_without_nc(self, var="biod"):
        if var == "biod":
            paths = self.test_algo(1)
        elif var == "fuel":
            paths = self.test_algo(2)
        screen1 = [x for x in list(paths.values()) if x < 1000]
        invert = [1 / (1 + x) for x in screen1]
        # replace 1/1001 or 1/101 by 0
        iic = sum(invert) / (self.size ** 4)
        return iic

    def CPL(self, var="biod"):
        if var == "biod":
            paths = self.test_algo(1)
        elif var == "fuel":
            paths = self.test_algo(2)
        # need to remove all pairs such that (i=j) from paths
        paths2 = [paths[k] for k in list(paths.keys()) if k[0] != k[1]]
        return stats.mean(paths2)

    def CPL_without_nc(self, var="biod"):
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
        if var == "zero":
            paths = self.test_algo(0)
        if var == "biod":
            paths = self.test_algo(1)
        elif var == "fuel":
            paths = self.test_algo(2)
        paths2 = [paths[x] for x in paths.keys() if paths[x] < 1000]
        return max(paths2)

    def landscape_shape_index(self, var="biod"):
        if var == 'biod':
            return 0.25*self.perimeter(var="biod")/self.size
        elif var == 'fuel':
            return 0.25*self.perimeter(var="fuel")/self.size

    def diagonal(self, direction, var="biod"):
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
        size = self.size - 1

        if var == "biod":
            threshold = 1
        elif var == "fuel":
            threshold = 2

        if i < 0 or j < 0 or i > size or j > size:
            return 0

        if i > 0:
            up = int(self.mat[i-1][j] == 0)
        elif i <= 0:
            up = 1

        if i < size:
            down = int(self.mat[i+1][j] == 0)
        elif i >= size:
            down = 1
        if j > 0:
            left = int(self.mat[i][j-1] == 0)
        elif j <= 0:
            left = 1
        if j < size:
            right = int(self.mat[i][j+1] == 0)
        elif j >= size:
            right = 1

        if i > 0 and j < size:
            up_right = int(self.mat[i-1][j+1] == 0)
        else:
            up_right = 1

        if i > 0 and j > 0:
            up_left = int(self.mat[i-1][j-1] == 0)
        else:
            up_left = 1

        if i < size and j > 0:
            down_left = int(self.mat[i+1][j-1] == 0)
        else:
            down_left = 1

        if i < size and j < size:
            down_right = int(self.mat[i+1][j+1] == 0)
        else:
            down_right = 1

        sum_neighbor = sum([up, down, right, left, up_right, up_left, down_right, down_left])
        return sum_neighbor*(self.mat[i][j] >= threshold)

    def nb_patches(self, var="biod"):
        size = self.size
        nb = []
        if var == "biod":
            for i in range(size):
                for j in range(size):
                    nb.append(self.number_zero_adjacent(i, j, var="biod"))
        elif var == "fuel":
            for i in range(size):
                for j in range(size):
                    nb.append(self.number_zero_adjacent(i, j, var="fuel"))

        return 1+sum([x == 8 for x in nb])
        # Components

    def components(self, var="biod"):
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
        if i < 0 or j < 0 or j > self.size - 1 or i > self.size - 1:
            return 0
        if var == "biod":
            if self.mat[i, j] == 0:
                return 0
            else:
                return 8 - self.number_zero_adjacent(i, j, var='biod')
        elif var == "fuel":
            if self.mat[i, j] == 0:
                return 0
            else:
                return 8 - self.number_zero_adjacent(i, j, var='fuel')

    def connectivity_correlation(self, var="biod", option_plot="No"):
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
        if var == 'biod':
            return self.perimeter(var="biod")/math.sqrt(sum(sum(self.mat == 1)))
        elif var == 'fuel':
            return self.perimeter(var="fuel")/math.sqrt(sum(sum(self.mat == 2)))

    def coord_zero(self):
        coord_zero = []
        for i in range(self.size):
            for j in range(self.size):
                if self.mat[i, j] == 0:
                    coord_zero.append((i, j))
        return coord_zero


def minDistance(grid, sourcex, sourcey, destx, desty):
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
    if ((x >= 0 and y >= 0) and
            (x < len(grid) and y < len(grid[0])) and
            (grid[x][y] != 0) and (visited[x][y] == False)):
        return True
    return False

#set of functions : for some reason, does not work with class
#shape = np.array((1,1,1,1,1,0,1,0,1,0,1,2,0,2,2,2))
def biod(shape):
    y = (shape >= 1).dot(utilities.connectivity_mat()).dot(np.transpose((shape >= 1)))
    return y
#biod(shape)

def biod_relative(shape):
    y = biod(shape)/biod(np.repeat(1,16))
    return y
#biod_relative(shape)

def fuel(shape):
    y = (shape >= 2).dot(utilities.connectivity_mat()).dot(np.transpose((shape >= 2)))
    return y
#fuel(shape)

def fuel_relative(shape):
    y = fuel(shape)/fuel(np.repeat(2,16))
    return y
#fuel_relative(shape)

def node_biod(shape):
    y = sum(shape>0)
    return y
#node_biod(shape)

def node_fuel(shape):
    y = sum(shape==2)
    return y
#node_fuel(shape)

def habitat_ratio(shape):
    y = sum(shape>0)/len(shape)
    return y
#habitat_ratio(shape)

def fuel_ratio(shape):
    y = node_fuel(shape)/len(shape)
    return y
#fuel_ratio(shape)

def test_algo(shape, k):
    if k >= 0:
        grid = utilities.to_matrix(shape>=k)
    elif k==0:
        grid = utilities.to_matrix(shape == k)
    grid = grid.astype(int)

    all_coordinates = []
    for i in range(int(math.sqrt(len(shape)))):
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
#convergence(succession)
def statistics_dataframe(data_path, output_=False):
    index = re.findall(r'\d+', data_path)
    #budget_here = index[2]
    path = "/home/simonjean/data/budget_" + str(params.budget) + "/successions_matched/"
    data = pd.read_csv(path+data_path).drop(columns=['Unnamed: 0'])

    data = data.iloc[:, 0:20]
    names = list(range(20))
    data.columns = names

    data["index"] = list(range(len(data)))
    # data['index_biod'] = [12]*len(data)
    # could be better to do in the end, with initial numbers as well.

    checker = data.melt(id_vars=["index"])
    checker["variable"] = checker['variable'].astype(int)
    checker = checker.sort_values(["index", "variable"])
    checker = checker.reset_index(drop=True)

    values = list(checker.value)
    values_unique = list(set(values))

    # region Setting up of data storage dictionnaries
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

    # region Compute statistics
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
        # print(values_unique.index(shape)/len(values_unique))

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

    checker.to_csv("/home/simonjean/data/budget_" + str(params.budget) + "/statistics/stats_"+data_path, index=False)
    if output_ == True:
        return checker.head()
def statistics_dataframe_light(data_path, output_=False):
    index = re.findall(r'\d+', data_path)
    #budget_here = index[2]
    path = "/home/simonjean/data/budget_" + str(params.budget) + "/successions_matched/"
    data = pd.read_csv(path+data_path).drop(columns=['Unnamed: 0'])

    data = data.iloc[:, 0:20]
    names = list(range(20))
    data.columns = names

    data["index"] = list(range(len(data)))
    # data['index_biod'] = [12]*len(data)
    # could be better to do in the end, with initial numbers as well.

    checker = data.melt(id_vars=["index"])
    checker["variable"] = checker['variable'].astype(int)
    checker = checker.sort_values(["index", "variable"])
    checker = checker.reset_index(drop=True)

    values = list(checker.value)
    values_unique = list(set(values))

    # region Setting up of data storage dictionnaries
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

    # region Compute statistics
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
        # print(values_unique.index(shape)/len(values_unique))

    checker["nodes_biod"] = [nodes_biod.get(n, n) for n in values]
    checker["nodes_fuel"] = [nodes_fuel.get(n, n) for n in values]
    checker["nodes_zero"] = [16 - checker.nodes_biod.loc[x] for x in range(len(checker))]
    checker["score_biod"] = [biod_score.get(n, n) for n in values]
    checker["score_fuel"] = [fuel_score.get(n, n) for n in values]

    checker["land_perimeter_biod"] = [land_perimeter_biod.get(n, n) for n in values]
    checker['land_perimeter_fuel'] = [land_perimeter_fuel.get(n, n) for n in values]
    checker["land_shape_index_biod"] = [land_shape_index_biod.get(n, n) for n in values]
    checker["land_shape_index_fuel"] = [land_shape_index_fuel.get(n, n) for n in values]

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

    checker.to_csv("/home/simonjean/data/budget_" + str(params.budget) + "/statistics/stats_"+data_path, index=False)
    if output_ == True:
        return checker.head()

def coord_potential_treatment(land):
    """
    Returns the coordinates of potential treatments for a landscape.
    :param land:
    :return:
    """
    counterfactual = np.minimum(land + np.repeat(1, params.size ** 2), np.repeat(2, params.size ** 2))
    coord_raw = np.where(counterfactual == 2)
    coord_raw = coord_raw[0].tolist()
    coord_treated = []
    for j in coord_raw:
        coord_treated.append((j // params.size, j % params.size))

    return coord_treated

def coord_treatment_2(land1,land2):
    '''
    Return where actual treatment was located in matrix coordinates. Associated with t, i.e, treatment applied in t leads to t+1
    :param land1: np.array
        land in t
    :param land2: np.array
        land in t+1
    :return:
    '''
    counterfactual = np.minimum(land1 + np.repeat(1,params.size**2), np.repeat(2, params.size**2))
    coord_raw = np.where(counterfactual != np.array(land2))
    coord_raw = coord_raw[0].tolist()
    coord_treated = []
    for j in coord_raw:
        coord_treated.append((j // params.size, j % params.size))

    return coord_treated

def distance_corner(coord_treated):
    """
    Define distance to corners of a grided landscape
    :param land1:
    :param land2:
    :param timing:
    :return:
    """
    dist = []
    for x in coord_treated:
        if x[0] <= 1 and x[1] <= 1:
            dist.append(math.sqrt(x[0] ** 2 + x[1] ** 2))
        elif x[0] <= 1 and x[1] > 1:
            dist.append(math.sqrt(x[0] ** 2 + ((params.size-1) - x[1]) ** 2))
        elif x[0] > 1 and x[1] <= 1:
            dist.append(math.sqrt(((params.size-1) - x[0]) ** 2 + x[1] ** 2))
        elif x[0] > 1 and x[1] > 1:
            dist.append(math.sqrt(((params.size-1) - x[0]) ** 2 + ((params.size-1) - x[1]) ** 2))

    return dist

def data_fix_coord(data_path):
    '''
    Fix data for coordinates of treatment. In row : potential treatment location, actual treatment employed to yield to next row
    :param data_path: int
        path to data
    :return: .csv
    '''
    data_dev = pd.read_csv("/home/simonjean/data/budget_" + str(params.budget) + "/statistics/" + data_path)

    #data_dev = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_"+str(params.budget)+"/"+data_path)

    data_process = data_dev["value"]
    data_process_l = list(data_process)
    data_process_l = [np.array(ast.literal_eval(x)) for x in data_process_l]

    potential_treatment = list(map(coord_potential_treatment, data_process_l))

    all_res = []
    lengths = []
    dist_realized = []
    dist_realized_mean = []
    dist_potential = []
    dist_potential_mean = []

    for x in range(0,len(data_process_l),20):
        a = []
        for row in range(19):
            a.append(coord_treatment_2(data_process_l[x+row],data_process_l[x+row+1]))
        a.append(None)
        all_res.extend(a)

    for y in all_res:
        try:
            lengths.append(len(y))
        except TypeError:
            lengths.append(None)

        if y==None:
            dist_realized.append(None)
            dist_realized_mean.append(None)
        else:
            try:
                dist_realized_mean.append(stats.mean(distance_corner(y)))
            except statistics.StatisticsError:
                dist_realized_mean.append(None)

            dist_realized.append(distance_corner(y))

    for z in potential_treatment:
        dist_potential.append(distance_corner(z))
        try:
            dist_potential_mean.append(stats.mean(distance_corner(z)))
        except statistics.StatisticsError:
            dist_potential_mean.append(None)

    data_dev['potential_treatment'] = potential_treatment
    data_dev['realized_treatment'] = all_res
    data_dev['nb_realized_treatment'] = lengths

    data_dev['average_distance_potential'] = dist_potential_mean
    data_dev['distance_potential'] = dist_potential
    data_dev['average_distance_realized'] = dist_realized_mean
    data_dev['distance_realized'] = dist_realized
    data_dev.to_csv("/home/simonjean/data/budget_" + str(params.budget) + "/statistics/" + data_path, index=False)

def candidates(opt='open'):
    candidates = os.listdir("/home/simonjean/data/budget_" + str(params.budget) + '/successions_matched')
    done = os.listdir("/home/simonjean/data/budget_" + str(params.budget) + "/statistics")

    done = [x.replace('stats_', '') for x in done]

    candidates_s = set(candidates)
    done_s = set(done)
    to_do = candidates_s - done_s
    candidates = list(to_do)

    if opt == 'write':
        file = open('/home/simonjean/data/budget_2/update.csv', 'w+', newline='')
    # writing the data into the file
        with file as f:
            write = csv.writer(f)
            write.writerows(candidates)

        return 'Candidates are written'

    elif opt == 'open':
        with open('/home/simonjean/data/budget_2/update.csv', 'r') as my_file:
            reader = csv.reader(my_file, delimiter = '\t')
            candidates = list(reader)

        candidates = [x[0].replace(',','') for x in candidates]
        return candidates


def data_fix(data_path):
    data = pd.read_csv('/home/simonjean/data/budget_' + str(params.budget) + "/statistics/" + data_path)

    values = list(data.value)
    values_unique = set(values)

    land_CPL_biod = {}
    land_CPL_fuel = {}
    land_shape_index_biod = {}
    land_shape_index_fuel = {}
    land_components_shape_index_biod = {}
    land_components_shape_index_fuel = {}
    for shape in values_unique:
        land = np.array(ast.literal_eval(shape))
        land_CPL_biod[shape] = CPL(land, nc=False)
        land_CPL_fuel[shape] = CPL(land, var="fuel", nc=False)
        land_shape_index_biod[shape] = landscape_shape_index(land)
        land_shape_index_fuel[shape] = landscape_shape_index(land, var="fuel")
        land_components_shape_index_biod[shape] = components_shape_index(land, var="biod")
        land_components_shape_index_fuel[shape] = components_shape_index(land, var="fuel")

    data["land_CPL_biod"] = [land_CPL_biod.get(n, n) for n in values]
    data["land_CPL_fuel"] = [land_CPL_fuel.get(n, n) for n in values]
    data["land_shape_index_biod"] = [land_shape_index_biod.get(n,n) for n in values]
    data['land_shape_index_fuel'] = [land_shape_index_fuel.get(n,n) for n in values]
    data['components_shape_index_biod'] = [land_components_shape_index_biod.get(n,n) for n in values]
    data['components_shape_index_fuel'] = [land_components_shape_index_fuel.get(n,n) for n in values]

    for x in range(20, len(data) + 1, 20):
        col = list(data.columns).index("value")
        succ = data.iloc[x - 20:x,col]
        conv = convergence(succ)
        data.iloc[x - 20:x, list(data.columns).index('convergence_time')] = conv[0]
        data.iloc[x - 20:x, list(data.columns).index('convergence_period')] = conv[1]

    data.to_csv('/home/simonjean/data/budget_' + str(params.budget) + "/statistics/" + data_path, index=False)

def enveloppe(data_path):
    data = pd.read_csv("/home/simonjean/data/budget_" + str(params.budget) + '/statistics/' + data_path)
    d_ = []
    data = data[['score_fuel','index_biod',"budget","convergence_time","variable"]]
    for x in range(0, len(data), 20):
        a = set(data.iloc[x:x + 20, data.columns.get_loc('convergence_time')]).pop()
        data_ = data.loc[x+a,:]
        data_ = list(data_[['score_fuel', 'index_biod', "budget"]])
        d_.append(data_)
    d_ = pd.DataFrame(d_, columns = ['score_fuel','index_biod',"budget"])
    d_ = d_.drop_duplicates()
    d_.to_csv("/home/simonjean/data/enveloppe/env_" + data_path, index=False)

def enveloppe2(data_path):
    data = pd.read_csv("/home/simonjean/data/budget_" + str(params.budget) + '/statistics/' + data_path)
    data = data[data['variable']>0]
    data = data[['score_fuel','index_biod',"budget"]]
    data = data.drop_duplicates()
    data.to_csv("/home/simonjean/data/enveloppe2/env_" + data_path, index=False)


def unique_convergence_patterns(data_path):
    '''
    Recover unique convergence patterns in whole data set
    :param data_path: str
        data path
    :return: .csv
    '''
    # Load data
    data = pd.read_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/' + data_path)
    # Set storage
    d_ = {}
    d_['patterns'] = set()
    # For every time series in the panel :
    for x in range(0, len(data), 20):
        # Find the convergence time and verify that it's unique
        a = set(data.iloc[x:x + 20, data.columns.get_loc('convergence_time')])
        if len(a) == 1:
            a = a.pop()
        else:
            print('Problemo dude')
        # Add to patterns only if not already in it (add and set())
        d_['patterns'].add(data.iloc[x + a, data.columns.get_loc('value')])
    d_['patterns'] = list(d_["patterns"])
    # Recover summary information
    d_['node_biod'] = [data.iloc[0, data.columns.get_loc('nodes_biod')]] * len(d_['patterns'])
    d_['node_fuel'] = [data.iloc[0, data.columns.get_loc('nodes_fuel')]] * len(d_['patterns'])
    d_['budget'] = [data.iloc[0, data.columns.get_loc('budget')]] * len(d_['patterns'])
    d_['constr_biod'] = [data.iloc[0, data.columns.get_loc('index_biod')]] * len(d_['patterns'])

    # Put dictionary in Dataframe format
    d_ = pd.DataFrame.from_dict(d_)
    # Save data
    d_.to_csv('/home/simonjean/data/patterns/' + data_path, index=False)

def matching_unique_landscape(data_path, keep="No"):
    """
    Match each landscape to its unique, non-equivalent counterpart.
    :param data_path: str
        data_path
    :param keep: str
        To store subsample, easier for download, exploratory analysis
    :return: .csv
    """
    # Load data
    data_check = pd.read_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/' + data_path)
    # Get rid of unnecessary data
    try:
        data_check = data_check.drop(['Unnamed:0','Unnamed:0.1'])
    except KeyError:
        pass
    # Load unique patterns data
    unique_patterns = pd.read_csv('/home/simonjean/data/unique_landscapes_convergence.csv')

    id_unique = []

    for x in range(0,int(len(data_check)),20):
        # For each landscape, find convergence time
        convergence_time = set(data_check.iloc[x:x+20,data_check.columns.get_loc('convergence_time')])
        if len(convergence_time)==1:
            convergence_time=convergence_time.pop()
        # Find unique pattern
        pattern = data_check.iloc[x+convergence_time, data_check.columns.get_loc('value')]
        # Match for unique  index (better for data storage volume)
        appendix = list(unique_patterns[unique_patterns['patterns']==pattern].iloc[0,:])[-1]
        id_unique.extend([appendix]*20)
    data_check['unique_id'] = id_unique
    data_check.to_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/' + data_path, index=False)
    # Second part :
    if keep == 'Yes':
        keep = pd.DataFrame()
        check = set(data_check.unique_id)
        check = list(check)
        count = []
        for x in check :
            count.append(sum((data_check[data_check['variable']==0].unique_id==x)))
        #Want to sample a little bit of the data to have a sense
        # Use index :
            sample_data = data_check[data_check['unique_id']==x]
            sample_from = list(sample_data['index'])
            to_sample = random.sample(sample_from, 4)
            new = data_check['index'].isin(to_sample)
            sample_data = data_check[new]
            keep = pd.concat([keep, sample_data])
        try :
            keep = keep.drop(columns=['Unnamed: 0', 'Unnamed: 0.1', 'Unnamed: 0.2'])
        except KeyError:
            pass
        keep.to_csv('/home/simonjean/data/budget_' + str(params.budget) + '/summary/sample/' + data_path, index=False)
    # We have kept some of the data


    # Now we need to do the procedure to keep the main results of the dataframe when it comes to the convergence landscape.
        # Keep data count
        # Should data count be kept ventilated by other variables?
            # keep indexes?

def re_indexing(datapath):
    """
    Fix to keep track of unique initial landscapes.
    :param datapath:
    :return:
    """
    files = os.listdir('/home/simonjean/data/all_data')
    d_ = pd.read_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/' + datapath)

    index = re.findall(r'\d+',datapath)
    index = index[2:]
    index.pop(-1)
    if len(index)==2:
        indexer = 'land4_' + index[0] + '_' + index[1] + '.pkl'
    elif len(index)==3:
        indexer = 'land4_' + index[0] + '_' + index[1] + 'cut_' + index[2] +'.pkl'

    index_d = list(d_['index'])
    index_d = [str(x) + '_' + str(files.index(indexer)) for x in index_d]
    d_['index_og'] = index_d
    d_.to_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/' + datapath, index=False)

def find_the_prob(data_path):
    a = [5, 26, 29, 43, 55, 63, 86, 95]
    d_ = pd.read_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/'+ data_path)
    saver = pd.DataFrame()
    for x in range(0,len(d_),20):
        if a.count(int(d_.iloc[x, d_.columns.get_loc('unique_id')]))>0:
            #print(int(d_.iloc[x, d_.columns.get_loc('unique_id')]))
            cv_time = set(d_.iloc[x:x+20,d_.columns.get_loc('convergence_time')]).pop()
            cv_period = set(d_.iloc[x:x+20, d_.columns.get_loc('convergence_period')]).pop()

            data_ = d_[(d_['variable'].isin(list(range(cv_time, cv_time+cv_period+1))) )& (d_['index']==int(x/20))]

            data_ = data_[['variable','value','convergence_time', 'convergence_period','unique_id']]
            saver = pd.concat([saver, data_])
    saver.to_csv('/home/simonjean/data/budget_' + str(params.budget) + '/troubleshooting/' + data_path)

def matching_cycles(datapath):
    '''
    Function for cycles attribution based on individual convergence landscape
    Uses cycles data
    :param datapath:
    :return:
    '''
    d_ = pd.read_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/' + datapath, engine='python')

    to_drop = d_.columns.to_list()
    to_drop = [x for x in to_drop if x.startswith('Unnamed')]
    d_ = d_.drop(columns=to_drop)

    cycle_index = []
    cycles = pd.read_csv('/home/simonjean/data/cycles.csv')

    for x in range(0, len(d_), 20):
        # Find convergence time and id
        convergence_time = set(d_.iloc[x:x + 20, d_.columns.get_loc('convergence_time')]).pop()
        id = set(d_.iloc[x:x + 20, d_.columns.get_loc('unique_id')]).pop()
        # Find what could be the successions :
        succession = set(cycles[cycles['phase1'] == id].phase2_pattern)
        succession = list(succession)
        succession = [tuple(ast.literal_eval(x)) for x in succession]
        # Find the index of potential successions
        index_c = set(cycles[cycles['phase1'] == id].cycle_index)
        index_c = list(index_c)
        # Find the actual data succession:
        candidate = d_.iloc[x + convergence_time + 1, d_.columns.get_loc('value')]
        candidate_formated = utilities.to_matrix(np.array(ast.literal_eval(candidate)))
        # For all potential successions following cycle initial phase
        for x in range(len(succession)):
            succession_s = set()
            #Define new set for potential succession
            succession_s.add(succession[x])
            # Check if all equivalents to actual succession might be equal to potential succession
            if len(set(utilities.equi_landscape_recover(candidate_formated)).intersection(succession_s)) == 1:
                # If so, assign cycle
                appendix = [index_c[x]] * 20
                cycle_index.extend(appendix)
    # Assign cycles to data
    d_['cycles'] = cycle_index
    # Save data
    d_.to_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/' + datapath, index=False)

def count_and_coordinates(data_path):
    """
    Count cycles and coordinates of t=0 and t=1 successions
    :param data_path:
    :return:
    """
    d_ = pd.read_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/' + data_path)

    # Data to keep
    budget = d_.iloc[0, d_.columns.get_loc('budget')]
    biod = d_.iloc[0, d_.columns.get_loc('index_biod')]
    node_biod_init = d_.iloc[0, d_.columns.get_loc('nodes_biod')]
    node_fuel_init = d_.iloc[0, d_.columns.get_loc('nodes_fuel')]

    # Count values
    a = d_['cycles'].value_counts()
    a = pd.concat([a,a])
    a = pd.DataFrame(a)
    a.reset_index(inplace=True)
    a = a.sort_values(by = ['index'])

    b = pd.DataFrame()
    index = []
    period = []
    for cycle in list(set(a['index'])):
        d__ = d_[d_['cycles']==cycle]
        v0 = sum([np.array(ast.literal_eval(x)) for x in d__.loc[d__['variable']==0, 'value']])
        v1 = sum([np.array(ast.literal_eval(x)) for x in d__.loc[d__['variable']==1, 'value']])
        index.extend([cycle]*2)
        period.extend([0,1])
        b = pd.concat([b, pd.DataFrame([np.transpose(v0), np.transpose(v1)])])
    del d__, d_

    b['index'] = index
    b['period'] = period

    a = pd.merge(a, b, on=['index'])
    a = a.drop_duplicates()

    # average
    a.loc[:,~a.columns.isin(['index', 'cycles', 'period'])] = a.loc[:, ~a.columns.isin(['index', 'cycles'])].div(a.cycles, axis=0)


    a['budget'] = [budget]*len(a)
    a['index_biod'] = [biod]*len(a)
    a['node_biod_init'] = [node_biod_init]*len(a)
    a['node_fuel_init'] = [node_fuel_init]*len(a)
    a = a[['node_biod_init','node_fuel_init', 'budget', 'index_biod', 'index','cycles', 0, 1, 2,
           3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]]
    a.to_csv('/home/simonjean/data/cycles_distrib/' + data_path, index=False)


def convergence_stats(datapath):
    d_ = pd.read_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/' + datapath)
    d_ = d_[['variable', 'nodes_biod', 'nodes_fuel', 'convergence_time', 'convergence_period', 'budget', 'index_biod',
             'score_fuel']]
    data = pd.DataFrame()
    for x in range(0, len(d_), 20):
        a = set(d_.iloc[x:x + 20, d_.columns.get_loc('convergence_time')]).pop()
        data_ = d_.loc[[x + a], :]
        data = pd.concat([data, data_], axis=0)

    data = data.drop_duplicates()
    data.to_csv('C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/convergences/' + datapath, index=False)

def cycle_constraint(datapath):
    data = pd.read_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/' + datapath)
    d_ = []
    index_biod = data.iloc[0,data.columns.get_loc('index_biod')]
    budget = data.iloc[0, data.columns.get_loc('budget')]
    data = data[[ "convergence_time","cycles"]]
    for x in range(0, len(data), 20):
        a = set(data.iloc[x:x + 20, data.columns.get_loc('convergence_time')]).pop()
        data_ = data.loc[x + a, :]
        data_ = list(data_[['cycles']])
        d_.append(data_)
    #data_ = data_.groupby(['index_biod', 'budget'])['cycles'].sum().reset_index()

    d_ = pd.DataFrame(d_, columns=['cycles'])
    d__ = pd.DataFrame(d_.cycles.value_counts()).reset_index()
    d__['index_biod'] = [index_biod]*len(d__)
    d__['budget'] = [budget]*len(d__)
    d__.to_csv('/home/simonjean/data/distribution_cycles/'+datapath)

def cv_time_period_by_cycle(datapath):
    d_ = pd.read_csv('/home/simonjean/data/budget_' + str(params.budget) + '/statistics/' + datapath)
    d_ = d_[[ 'convergence_time', 'convergence_period', 'budget', 'index_biod','cycles']]
    data = pd.DataFrame()
    for x in range(0, len(d_), 20):
        a = set(d_.iloc[x:x + 20, d_.columns.get_loc('convergence_time')]).pop()
        data_ = d_.loc[[x + a], :]
        data = pd.concat([data, data_], axis=0)

    data = data.groupby(list(data.columns)).size().reset_index(name='count')
    data = data.drop_duplicates()
    data.to_csv('/home/simonjean/data/convergences/convergence_april/' + datapath, index=False)
