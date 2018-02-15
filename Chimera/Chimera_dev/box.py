import sys
import pandas as pd
import numpy as np
import time
from Chimera.Chimera_dev import mesh, console, backends, neighbors_cy
# from . import neighbors
import warnings
warnings.filterwarnings('ignore')
import pyximport; pyximport.install()


class Box:

    def __init__(self, evolution_time, conduction=True, settling_mode='stokes terminal',
                 radioactivity=True, chemistry=True, verbose=True):
        self.mesh = pd.DataFrame({
        })
        self.objects = pd.DataFrame({
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


    def box_info(self):
        df_memory = self.mesh.memory_usage(deep=True)
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
        self.spatial_sigfigs = m.get_spatial_sigfigs()
        console.nominal("Assigning mesh to dataframe...", verbose=self.verbose)
        m.build(df=self.mesh)
        console.nominal("Exiting box construction...", verbose=self.verbose)
        console.event("Box constructed!", verbose=self.verbose)
        console.event("Fetching nearest neighbors...", verbose=self.verbose)
        self.mesh['object'] = ['NONE' for i in range(len(self.mesh['coords']))]
        self.mesh['object_id'] = np.NAN
        self.mesh['neighbors'] = np.NAN
        self.mesh['xplus_index'] = [0 for i in range(len(self.mesh['coords']))]
        self.mesh['xminus_index'] = [0 for i in range(len(self.mesh['coords']))]
        self.mesh['yplus_index'] = [0 for i in range(len(self.mesh['coords']))]
        self.mesh['yminus_index'] = [0 for i in range(len(self.mesh['coords']))]
        self.mesh['zplus_index'] = [0 for i in range(len(self.mesh['coords']))]
        self.mesh['zminus_index'] = [0 for i in range(len(self.mesh['coords']))]
        self.mesh['temperature'] = [0.0 for i in range(len(self.mesh['coords']))]
        self.objects['object'] = np.NAN
        self.objects['object_id'] = np.NAN
        self.objects['curr_loc'] = np.NAN
        self.objects['temperature'] = np.NAN
        self.objects['radius'] = np.NAN
        spatial_sigfigs = m.get_spatial_sigfigs()
        neighbor_count = 0
        total_count = len(self.mesh['coords'])
        t_start = time.time()
        arr = self.mesh['coords'].tolist()
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
        self.mesh['xplus_index'] = x_plus
        self.mesh['xminus_index'] = x_minus
        self.mesh['yplus_index'] = y_plus
        self.mesh['yminus_index'] = y_minus
        if self.dimension is 3:
            self.mesh['zplus_index'] = z_plus
            self.mesh['zminus_index'] = z_minus
        console.event("Finished finding nearest neighbors! (task took {}s for {} points)".format(
                        time.time() - t_start, total_count),
                      verbose=self.verbose)
        return self.mesh

    def insert_matrix(self, material, temperature, depth_range):
        console.event("Inserting matrix ({}) into the box!".format(material), verbose=self.verbose)
        t_start = time.time()
        subdf = self.mesh.query('{} <= z_vals <= {}'.format(depth_range[0], depth_range[1]))
        self.mesh.drop(['z_vals'], axis=1)
        coords = subdf['coords'].tolist()
        temperatures = self.mesh['temperature'].tolist()
        objects = self.mesh['object'].tolist()
        object_ids = self.mesh['object_id'].tolist()
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
        self.mesh['temperature'] = temperatures
        self.mesh['object'] = objects
        self.mesh['object_id'] = object_ids
        self.matrix = True
        console.event("Finished inserting matrix ({}) into the box! (task took {}s)".format(material,
                                                                                            time.time() - t_start), verbose=self.verbose)

    def insert_boundary(self, temperature, depth_range, material='boundary'):
        console.event("Inserting boundary ({}) into the box!".format(material), verbose=self.verbose)
        t_start = time.time()
        subdf = self.mesh.query('{} <= z_vals <= {}'.format(depth_range[0], depth_range[1]))
        coords = subdf['coords'].tolist()
        temperatures = self.mesh['temperature'].tolist()
        objects = self.mesh['object'].tolist()
        object_ids = self.mesh['object_id'].tolist()
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
        self.mesh['temperature'] = temperatures
        self.mesh['object'] = objects
        self.mesh['object_id'] = object_ids
        self.boundary = True
        console.event("Finished inserting boundary ({}) into the box! (task took {}s)".format(material,
                                                                                              time.time() - t_start), verbose=self.verbose)


    def insert_object(self, material, temperature, radius, x, y, z=None):
        console.event("Inserting object ({}) into the box!".format(material), verbose=self.verbose)
        if z is None and self.dimension is 3:
            console.error("Box dimension is 3 and z is not specified!", verbose=self.verbose)
            sys.exit(1)
        if self.dimension is 3 and z is not None:
            requested_coord = (x, y, z)
        elif self.dimension is 2 and z is None:
            requested_coord = (x, y)
        t_start = time.time()
        objects = self.objects['object'].tolist()
        object_ids = self.objects['object_id'].tolist()
        temperatures = self.objects['temperature'].tolist()
        radii = self.objects['radius'].tolist()
        locs = self.objects['curr_loc'].tolist()
        console.nominal("Inserting object ({}) at {}...".format(material, requested_coord), verbose=self.verbose)
        objects.append(material)
        object_ids.append(backends.generate_object_id(object_type='object',
                                                      id_val=self.id_val))
        radii.append(radius)
        locs.append(requested_coord)
        temperatures.append(temperature)
        self.objects['object'] = objects
        self.objects['object_id'] = object_ids
        self.objects['temperature'] = temperatures
        self.objects['radius'] = radii
        self.objects['curr_loc'] = locs
        self.id_val += 1
        self.object = True
        console.event("Finished inserting object ({}) into the box! (task took {}s)".format(material,
                                                                                            time.time() - t_start), verbose=self.verbose)


    def verify_box(self):
        objects = self.mesh['object'].tolist()
        if 'NONE' in objects:
            console.error("WARNING! Not all coordinate spaces defined in box!", verbose=True)
            self.mesh.drop(['z_vals'])
        else:
            console.event("Box verified!", verbose=True)


    def to_csv(self):
        self.mesh.to_csv("mesh.csv")
        self.objects.to_csv("objects.csv")


    def update(self, time=1, conduction=True, settling='stokes terminal'):
        pass
        
