import os
import sys
import pandas as pd
import numpy as np
import time
from . import console
from . import mesh
from . import backends
# from . import neighbors
import warnings
warnings.filterwarnings('ignore')
import pyximport; pyximport.install()
from . import neighbors_cy

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
        self.dimension = None
        self.matrix = None
        self.object = None
        self.boundary = None
        self.id_val = 0

    def update(self):
        pass

    def box_info(self):
        df_memory = self.space.memory_usage(deep=True)
        return df_memory

    def build(self, spatial_res, x, y, z=None):
        if z is not 3:
            self.dimension = 3
            self.max_x, self.max_y, self.max_z = x, y, z
        else:
            self.dimension = 2
            self.max_x, self.max_y, self.max_z = x, y, None
        self.spatial_res = spatial_res
        console.event("Constructing box...", verbose=self.verbose)
        console.nominal("Building mesh...", verbose=self.verbose)
        m = mesh.Mesh(spatial_res=spatial_res, x=x, y=y, z=z, verbose=self.verbose)
        console.nominal("Assigning mesh to dataframe...", verbose=self.verbose)
        m.build(df=self.space)
        console.nominal("Exiting box construction...", verbose=self.verbose)
        console.event("Box constructed!", verbose=self.verbose)
        console.event("Fetching nearest neighbors...", verbose=self.verbose)
        self.space['object'] = ['NONE' for i in range(len(self.space['coords']))]
        self.space['object_id'] = np.NAN
        self.space['neighbors'] = np.NAN
        self.space['xplus_index'] = [0 for i in range(len(self.space['coords']))]
        self.space['xminus_index'] = [0 for i in range(len(self.space['coords']))]
        self.space['yplus_index'] = [0 for i in range(len(self.space['coords']))]
        self.space['yminus_index'] = [0 for i in range(len(self.space['coords']))]
        self.space['zplus_index'] = [0 for i in range(len(self.space['coords']))]
        self.space['zminus_index'] = [0 for i in range(len(self.space['coords']))]
        self.space['temperature'] = [0 for i in range(len(self.space['coords']))]
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
            if self.dimension is 3:
                console.nominal("Finding neighbor for ({}, {}, {}) ({}/{} points)".format(
                    coords[0], coords[1], coords[2], neighbor_count, total_count), verbose=self.verbose)
            else:
                console.nominal("Finding neighbor for ({}, {})".format(
                    coords[0], coords[1]), verbose=self.verbose)
            n = neighbors_cy.get_neighbors(verbose=self.verbose, coords=coords, array=arr, max_x=x, max_y=y, max_z=z,
                                               spatial_res=spatial_res, spatial_sigfigs=spatial_sigfigs)
            x_plus.append(n[0])
            x_minus.append(n[1])
            y_plus.append(n[2])
            y_minus.append(n[3])
            if self.dimension is 3:
                z_plus.append(n[4])
                z_minus.append(n[5])
            neighbor_count += 1
        self.space['xplus_index'] = x_plus
        self.space['xminus_index'] = x_minus
        self.space['yplus_index'] = y_plus
        self.space['yminus_index'] = y_minus
        if self.dimension is 3:
            self.space['zplus_index'] = z_plus
            self.space['zminus_index'] = z_minus
        console.event("Finished finding nearest neighbors! (task took {}s for {} points)".format(
                        time.time() - t_start, total_count),
                      verbose=self.verbose)
        return self.space

    def insert_matrix(self, material, temperature, depth_range):
        console.event("Inserting matrix ({}) into the box!".format(material), verbose=self.verbose)
        t_start = time.time()
        subdf = self.space.query('{} <= z_vals <= {}'.format(depth_range[0], depth_range[1]))
        self.space.drop(['z_vals'], axis=1)
        coords = subdf['coords'].tolist()
        temperatures = self.space['temperature'].tolist()
        objects = self.space['object'].tolist()
        object_ids = self.space['object_id'].tolist()
        len_coords = len(coords)
        for coord in coords:
            if depth_range[0] <= coords[2][self.dimension - 1] <= depth_range[1]:
                index = backends.predict_index(coord=coord, max_x=self.max_x, max_y=self.max_y, max_z=self.max_z,
                                               spatial_res=self.spatial_res, verbose=self.verbose)
                console.nominal("Inserting matrix ({}) at {}...".format(material, coord), verbose=self.verbose)
                temperatures[index] = temperature
                objects[index] = material
                object_ids[index] = backends.generate_object_id(object_type='matrix',
                                                                id_val=self.id_val)
                self.id_val += 1
        self.space['temperature'] = temperatures
        self.space['object'] = objects
        self.space['object_id'] = object_ids
        self.matrix = True
        console.event("Finished inserting matrix ({}) into the box! (task took {}s)".format(material,
                        time.time() - t_start), verbose=self.verbose)

    def insert_boundary(self, temperature, depth_range, material='boundary'):
        console.event("Inserting boundary ({}) into the box!".format(material), verbose=self.verbose)
        t_start = time.time()
        subdf = self.space.query('{} <= z_vals <= {}'.format(depth_range[0], depth_range[1]))
        coords = subdf['coords'].tolist()
        temperatures = self.space['temperature'].tolist()
        objects = self.space['object'].tolist()
        object_ids = self.space['object_id'].tolist()
        len_coords = len(coords)
        for coord in coords:
            if depth_range[0] <= coords[2][self.dimension - 1] <= depth_range[1]:
                index = backends.predict_index(coord=coord, max_x=self.max_x, max_y=self.max_y, max_z=self.max_z,
                                               spatial_res=self.spatial_res, verbose=self.verbose)
                console.nominal("Inserting boundary ({}) at {}...".format(material, coord), verbose=self.verbose)
                temperatures[index] = temperature
                objects[index] = material
                object_ids[index] = backends.generate_object_id(object_type='boundary',
                                                                id_val=self.id_val)
                self.id_val += 1
        self.space['temperature'] = temperatures
        self.space['object'] = objects
        self.space['object_id'] = object_ids
        self.boundary = True
        console.event("Finished inserting boundary ({}) into the box! (task took {}s)".format(material,
                        time.time() - t_start), verbose=self.verbose)


    def insert_object(self, material, temperature, x, y, z=None):
        console.event("Inserting object ({}) into the box!".format(material), verbose=self.verbose)
        if z is None and self.dimension is 3:
            console.error("Box dimension is 3 and z is not specified!", verbose=self.verbose)
            sys.exit(1)
        if self.dimension is 3 and z is not None:
            requested_coord = (x, y, z)
        elif self.dimension is 2 and z is None:
            requested_coord = (x, y)
        t_start = time.time()
        coords = self.space['coords'].tolist()
        temperatures = self.space['temperature'].tolist()
        objects = self.space['object'].tolist()
        object_ids = self.space['object_id'].tolist()
        len_coords = len(coords)
        index = coords.index(requested_coord)
        console.nominal("Inserting object ({}) at {}...".format(material, requested_coord), verbose=self.verbose)
        temperatures[index] = temperature
        objects[index] = material
        object_ids[index] = backends.generate_object_id(object_type='object',
                                                                id_val=self.id_val)
        self.id_val += 1
        self.space['temperature'] = temperatures
        self.space['object'] = objects
        self.space['object_id'] = object_ids
        self.object = True
        console.event("Finished inserting object ({}) into the box! (task took {}s)".format(material,
                      time.time() - t_start), verbose=self.verbose)


    def verify_box(self):
        objects = self.space['object'].tolist()
        if 'NONE' in objects:
            console.error("WARNING! Not all coordinate spaces defined in box!", verbose=True)
        else:
            console.event("Box verified!", verbose=True)



    def to_csv(self):
        self.space.to_csv("space.csv")
