import os
import numpy as np
import pandas as pd
import time
from . import console


class Mesh:

    def __init__(self, spatial_res, x, y, z=None, verbose=True):
        self.verbose = verbose
        self.spatial_res = spatial_res
        if '.' in str(self.spatial_res):
            self.spatial_sigfigs = len(str(self.spatial_res)) - 1
        else:
            self.spatial_sigfigs = len(str(self.spatial_res))
        self.x = x
        self.y = y
        self.x_coords = None
        self.y_coords = None
        if z is not None:
            self.z = z
            self.z_coords = None
            self.dimensions = 3
        else:
            self.dimensions = 2



    def build(self, df=None):
        console.event("Generating model vertices...", verbose=self.verbose)
        t_start = time.time()
        x_coords = np.arange(0, self.x + self.spatial_res, self.spatial_res)
        y_coords = np.arange(0, self.y + self.spatial_res, self.spatial_res)
        if self.dimensions is 3:
            console.nominal("Detected a 3-D system! Generating a {} x {} x {} model!".format(
                round(self.x / self.spatial_res, self.spatial_sigfigs),
                round(self.y / self.spatial_res, self.spatial_sigfigs),
                round(self.z / self.spatial_res, self.spatial_sigfigs)),
                verbose=self.verbose)
            z_coords = np.arange(0, self.z + self.spatial_res, self.spatial_res)
            mesh = np.meshgrid(x_coords, y_coords, z_coords)
            nodes = list(zip(*(dim.flat for dim in mesh)))
            cleaned_nodes = []
            for i in nodes:
                x = round(i[0], self.spatial_sigfigs)
                y = round(i[1], self.spatial_sigfigs)
                z = round(i[2], self.spatial_sigfigs)
                cleaned_nodes.append((x, y, z))
        else:
            console.nominal("Detected a 2-D system! Generating a {} x {} model!".format(
                round(self.x / self.spatial_res, self.spatial_sigfigs),
                round(self.y / self.spatial_res, self.spatial_sigfigs)),
                verbose=self.verbose)
            mesh = np.meshgrid(x_coords, y_coords)
            nodes = list(zip(*(dim.flat for dim in mesh)))
            cleaned_nodes = []
            for i in nodes:
                x = round(i[0], self.spatial_sigfigs)
                y = round(i[1], self.spatial_sigfigs)
                cleaned_nodes.append((x, y))
        if df is not None:
            df['coords'] = (str(list(i)) for i in nodes)
        time_task = time.time() - t_start
        console.event("Finished generating model vertices! (task took {}s)".format(time_task), verbose=self.verbose)
        return cleaned_nodes
