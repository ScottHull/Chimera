import os
import pandas as pd
import numpy as np
import time
from . import console
from . import mesh
from . import neighbors
import warnings
warnings.filterwarnings('ignore')

class Box:

    def __init__(self, evolution_time, conduction=True, settling_mode='stokes terminal',
                 radioactivity=True, chemistry=True, verbose=True):
        self.space = pd.DataFrame({
        })
        self.conduction = conduction
        self.settling_mode = settling_mode
        self.radioactivity = radioactivity
        self.chemistry = chemistry
        self.initial_time = evolution_time
        self.evolution_time = evolution_time
        self.verbose = verbose
        self.max_x = 0.0
        self.max_y = 0.0
        self.max_z = 0.0
        self.spatial_res = 0.0
        self.spatial_sigfigs = 0.0

    def update(self):
        pass

    def box_info(self):
        df_memory = self.space.memory_usage(deep=True)
        return df_memory

    def build(self, spatial_res, x, y, z=None):
        console.event("Constructing box...", verbose=self.verbose)
        console.nominal("Building mesh...", verbose=self.verbose)
        m = mesh.Mesh(spatial_res=spatial_res, x=x, y=y, z=z, verbose=self.verbose)
        console.nominal("Assigning mesh to dataframe...", verbose=self.verbose)
        m.build(df=self.space)
        console.nominal("Exiting box construction...", verbose=self.verbose)
        console.event("Box constructed!", verbose=self.verbose)
        console.event("Fetching nearest neighbors...", verbose=self.verbose)
        self.space['neighbors'] = np.NAN
        self.space['xplus_index'] = [0 for i in range(len(self.space['coords']))]
        self.space['xminus_index'] = [0 for i in range(len(self.space['coords']))]
        self.space['yplus_index'] = [0 for i in range(len(self.space['coords']))]
        self.space['yminus_index'] = [0 for i in range(len(self.space['coords']))]
        self.space['zplus_index'] = [0 for i in range(len(self.space['coords']))]
        self.space['zminus_index'] = [0 for i in range(len(self.space['coords']))]
        spatial_sigfigs = m.get_spatial_sigfigs()
        neighbor_count = 1
        total_count = len(self.space['coords'])
        t_start = time.time()
        arr = self.space['coords'].tolist()
        x_plus = []
        x_minus = []
        y_plus = []
        y_minus = []
        z_plus = []
        z_minus = []
        for coords in arr:
            if z is not None:
                console.nominal("Finding neighbor for ({}, {}, {}) ({}/{} points)".format(
                    coords[0], coords[1], coords[2], neighbor_count, total_count), verbose=self.verbose)
            else:
                console.nominal("Finding neighbor for ({}, {})".format(
                    coords[0], coords[1]), verbose=self.verbose)
            n = neighbors.get_neighbors(verbose=self.verbose, coords=coords, array=arr, max_x=x, max_y=y, max_z=z,
                                               spatial_res=spatial_res, spatial_sigfigs=spatial_sigfigs)
            x_plus.append(n[0])
            x_minus.append(n[1])
            y_plus.append(n[2])
            y_minus.append(n[3])
            if z is not None:
                z_plus.append(n[4])
                z_minus.append(n[5])
            neighbor_count += 1
        self.space['xplus_index'] = x_plus
        self.space['xminus_index'] = x_minus
        self.space['yplus_index'] = y_plus
        self.space['yminus_index'] = y_minus
        if z is not None:
            self.space['zplus_index'] = z_plus
            self.space['zminus_index'] = z_minus
        console.event("Finished finding nearest neighbors! (task took {}s for {} points)".format(
                        time.time() - t_start, total_count),
                      verbose=self.verbose)
        return self.space