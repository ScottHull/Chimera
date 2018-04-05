import os
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors
import matplotlib.colorbar
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from itertools import product, combinations
import moviepy.editor as mpy
import numpy as np
import time
import shutil
from . import console

class plots():

    def __init__(self):
        self.plot_cell_frames = []
        self.plot_cell_save = False
        self.temperatures = None
        self.colorsmap = None
        self.norm_colors = None
        if "plot_cell" in os.listdir(os.getcwd()):
            shutil.rmtree(os.getcwd() + "/plot_cell")
        os.mkdir(os.getcwd() + "/plot_cell")

    def plot_cell(self, object_coords, nearest_coords, vertex_indices, mesh_coords, max_x, max_y, max_z, spatial_res,
                  model_time, temperatures, heat=False, save=False, show=False):
        if save is True or show is True:
            fig = plt.figure()
            ax = Axes3D(fig)
            if heat is True:
                if self.temperatures is None:
                    self.temperatures = temperatures
                    self.norm_colors = mpl.colors.Normalize(vmin=2000, vmax=3500)
                    self.colorsmap = matplotlib.cm.ScalarMappable(norm=self.norm_colors, cmap='jet')
                    self.colorsmap.set_array(temperatures)
                # norm_colors = mpl.colors.Normalize(vmin=min(temperatures), vmax=max(temperatures))
                # norm_colors = mpl.colors.Normalize(vmin=2000, vmax=3500)
                # colorsmap = matplotlib.cm.ScalarMappable(norm=norm_colors, cmap='jet')
                # colorsmap.set_array(temperatures)
                # cb = fig.colorbar(colorsmap)
                cb = fig.colorbar(self.colorsmap)
                ax.scatter(*zip(*mesh_coords), marker='s', s=5, c=temperatures, cmap='jet', alpha=0.25)
            for index, object_coord in enumerate(object_coords):
                x, y, z = object_coord[0], object_coord[1], object_coord[2]
                cell_vertices = []
                for i in vertex_indices[index]:
                    cell_vertices.append(list(mesh_coords[i]))
                ax.scatter3D(x, y, z, color='r')
                points = np.array(cell_vertices)
                ax.scatter(points[:, 0], points[:, 1], points[:, 2], alpha=0.5, s=0.5)
                for s, e in combinations(points, 2):
                    # All diagonals will be greater than 2
                    if np.sum(np.abs(s - e)) <= spatial_res:
                        ax.plot3D(*zip(s, e), color="k", alpha=0.5)
                # nearest_coord = nearest_coords[index]
                # ax.scatter3D(mesh_coords[nearest_coord][0], mesh_coords[nearest_coord][1], mesh_coords[nearest_coord][2], color='y')
            ax.set_title("Model time: {}".format(model_time))
            ax.set_xlim(xmin=0.0, xmax=max_x)
            ax.set_ylim(ymin=0.0, ymax=max_y)
            ax.set_zlim(zmin=0.0, zmax=max_z)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            ax.invert_zaxis()
            if save is True:
                fig.savefig(os.getcwd() + "/plot_cell/" + str(model_time) + ".png", formmat='png')
                self.plot_cell_frames.append(str(model_time) + ".png")
                self.plot_cell_save = True
            if show is True:
                plt.show()

    def animate(self, initial_time):
        if self.plot_cell_save is True:
            home_dir = os.getcwd()
            frames = self.plot_cell_frames
            dir = os.getcwd() + "/plot_cell"
            os.chdir(dir)
            animation = mpy.ImageSequenceClip(frames, fps=(initial_time / (initial_time / 3)), load_images=True)
            os.chdir(home_dir)
            animation.write_gif('plot_cell.gif', fps=(initial_time / (initial_time / 3)))