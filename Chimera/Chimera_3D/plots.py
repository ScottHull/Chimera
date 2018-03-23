import os
import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import numpy as np
import time
from . import console

def plot_cell(object_coord, nearest_coord, vertex_indeces, mesh_coords):
    # print(vertex_indeces)
    fig = plt.figure()
    ax = Axes3D(fig)
    x, y, z = object_coord[0], object_coord[1], object_coord[2]
    cell_verteces = []
    for i in vertex_indeces:
        cell_verteces.append(list(mesh_coords[i]))
    ax.scatter3D(x, y, z, color='r')
    for i in cell_verteces:
        x, y, z = i[0], i[1], i[2]
        ax.scatter3D(x, y, z, color='b')
    ax.scatter3D(mesh_coords[nearest_coord][0], mesh_coords[nearest_coord][1], mesh_coords[nearest_coord][2], color='y')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.show()
