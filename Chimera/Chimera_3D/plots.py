import os
import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import numpy as np
import time
from . import console

# def plot_cell(object_coord, vertex_indeces, mesh_coords):
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     r = [-1, 1]
#     cell_verteces = []
#     for i in vertex_indeces:
#         cell_verteces.append(list(mesh_coords[i]))
#     points = np.array(cell_verteces)
#     X, Y = np.meshgrid(r, r)
#     ax.plot_surface(X, Y, 1, alpha=0.5)
#     ax.plot_surface(X, Y, -1, alpha=0.5)
#     ax.plot_surface(X, -1, Y, alpha=0.5)
#     ax.plot_surface(X, 1, Y, alpha=0.5)
#     ax.plot_surface(1, X, Y, alpha=0.5)
#     ax.plot_surface(-1, X, Y, alpha=0.5)
#     ax.scatter3D(points[:, 0], points[:, 1], points[:, 2])
#     fig.show()