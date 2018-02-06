import os
import numpy as np
import pandas as pd


class Mesh:

    def __init__(self, spatial_res, x, y, z=None):
        self.space = pd.DataFrame({
        })
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



    def build(self, apply_to_df=True):
        coords = []
        self.x_coords = np.arange(0, self.x, round(self.spatial_res, self.spatial_sigfigs))
        self.y_coords = np.arange(0, self.y, round(self.spatial_res, self.spatial_sigfigs))
        if self.dimensions is 3:
            self.z_coords = np.arange(0, self.z, round(self.spatial_res, self.spatial_sigfigs))
            for x in self.x_coords:
                for y in self.y_coords:
                    for z in self.z_coords:
                        coord_set = (x, y, z)
                        coords.append(coord_set)
        else:
            for x in self.x_coords:
                for y in self.y_coords:
                    coord_set = (x, y)
                    coords.append(coord_set)
        if apply_to_df is True:
            self.space['coords'] = (i for i in coords)
