from . import mesh
from . import console
from . import neighbors
from . import backends
from . import heat
import os
import pandas as pd
import numpy as np
import time

class Line:
    """
    Handles and tracks all functions with the construction and evolution of the linear Chimera model.
    """
    def __init__(self, evolution_time, time_units='s', conduction=True, convection=True, settling_mode='stokes terminal',
                 radioactivity=True, chemistry=True, verbose=True):
        self.mesh = pd.DataFrame({
            'coords': [],
            'object': [],
            'object_id': [],
            'zplus_index': [],
            'zminus_index': [],
            'temperature': [],
            'conductivity': [],
            'dT_dt': [],
            'density': [],
            'viscosity': []
        })
        self.objects = pd.DataFrame({
            'object': [],
            'object_id': [],
            'coords': [],
            'density': [],
            'radius': [],
            'velocity': [],
        })
        self.layers = pd.DataFrame({
            'object': [],
            'min_z': [],
            'max_z': [],
            'min_z_temp': [],
            'max_z_temp': [],
            'min_z_dT_dt': [],
            'max_z_dT_dt': [],
            'nusselt': [],
        })
        self.conduction = conduction
        self.convection = convection
        self.settling_mode = settling_mode
        self.radioactivity = radioactivity
        self.chemistry = chemistry
        self.initial_time = evolution_time
        self.evolution_time = evolution_time
        self.verbose = verbose
        self.max_z = 0.0
        self.spatial_res = 0.0
        self.spatial_sigfigs = 0.0
        self.delta_time = None
        self.time_units = None
        if time_units[0].lower() == 's':
            self.time_units = 'seconds'
        else:
            self.time_units = 'year'
        self.dimension = None
        self.matrix = None
        self.object = None
        self.boundary = None
        self.id_val = 0
        self.materials = {}
        self.conductivities = []

    def build(self, spatial_res, z):
        """
        Construct the mesh of the model and instantiate the mesh dataframe.
        :param spatial_res:
        :param z:
        :return:
        """
        self.spatial_res = spatial_res
        self.max_z = z
        console.event("Constructing line...", verbose=self.verbose)
        console.nominal("Building mesh...", verbose=self.verbose)
        m = mesh.Mesh(spatial_res=spatial_res, z=z, verbose=self.verbose)
        self.spatial_sigfigs = m.get_spatial_sigfigs()
        console.nominal("Assigning mesh to dataframe...", verbose=self.verbose)
        m.build_linear(df=self.mesh)
        console.nominal("Exiting line construction...", verbose=self.verbose)
        console.event("Line constructed!", verbose=self.verbose)
        self.mesh['object'] = ['NONE' for i in range(len(self.mesh['coords']))]
        self.mesh['object_id'] = np.NAN
        self.mesh['zplus_index'] = [0 for i in range(len(self.mesh['coords']))]
        self.mesh['zminus_index'] = [0 for i in range(len(self.mesh['coords']))]
        self.mesh['temperature'] = [0.0 for i in range(len(self.mesh['coords']))]
        self.mesh['conductivity'] = [0.0 for i in range(len(self.mesh['coords']))]
        self.mesh['dT_dt'] = [0.0 for i in range(len(self.mesh['coords']))]
        self.mesh['density'] = [0.0 for i in range(len(self.mesh['coords']))]
        self.mesh['viscosity'] = [0.0 for i in range(len(self.mesh['coords']))]
        console.event("Fetching nearest neighbors...", verbose=self.verbose)
        neighbor_count = 1
        total_count = len(self.mesh['coords'])
        t_start = time.time()
        arr = self.mesh['coords'].tolist()
        z_plus = []
        z_minus = []
        for coords in arr:
            console.nominal("Finding neighbor for ({}) ({}/{} points)".format(
                coords, neighbor_count, total_count), verbose=self.verbose)
            n = neighbors.get_linear_neighbors(verbose=self.verbose, coord=coords, max_z=self.max_z,
                                               spatial_sigfigs=self.spatial_sigfigs, spatial_res=self.spatial_res)
            z_plus.append(n[0])
            z_minus.append(n[1])
            neighbor_count += 1
        self.mesh['zplus_index'] = z_plus
        self.mesh['zminus_index'] = z_minus
        console.event("Finished finding nearest neighbors! (task took {}s for {} points)".format(
            time.time() - t_start, total_count),
            verbose=self.verbose)
        return self.mesh

    def insert_matrix(self, material, conductivity, initial_temp, initial_density,
                      viscosity, depth_range, temp_grad=None):
        """
        Insert a matrix material into the linear model within a specified depth range.
        :param material:
        :param conductivity:
        :param initial_temp:
        :param depth_range:
        :param temp_grad:
        :return:
        """
        console.event("Inserting matrix ({}) into the model!".format(material), verbose=self.verbose)
        t_start = time.time()
        subdf = self.mesh.query('{} <= coords <= {}'.format(depth_range[0], depth_range[1]))
        coords = subdf['coords'].tolist()
        temperatures = self.mesh['temperature'].tolist()
        layer_temperatures = []
        objects = self.mesh['object'].tolist()
        object_ids = self.mesh['object_id'].tolist()
        conductivities = self.mesh['conductivity'].tolist()
        temperature = initial_temp
        for coord in coords:
            if depth_range[0] <= coord <= depth_range[1]:
                index = backends.predict_linear_index(z=coord, spatial_res=self.spatial_res)
                console.nominal("Inserting matrix ({}) at {}...".format(material, coord), verbose=self.verbose)
                temperatures[index] = temperature
                layer_temperatures.append(temperature)
                conductivities[index] = conductivity
                objects[index] = material
                object_ids[index] = backends.generate_object_id(object_type='matrix',
                                                                id_val=self.id_val)
                self.id_val += 1
                if temp_grad is not None:
                    temperature += temp_grad
        self.mesh['temperature'] = temperatures
        self.mesh['object'] = objects
        self.mesh['object_id'] = object_ids
        self.mesh['conductivity'] = conductivities
        if material not in self.materials.keys():
            self.materials.update({material: float(conductivity)})
            self.conductivities.append(conductivity)
        layers_row = {'object': material,
                      'min_z': depth_range[0],
                      'max_z': depth_range[1],
                      'z_min_temp': layer_temperatures[0],
                      'z_max_temp': layer_temperatures[-1],
                      'min_z_dT_dt': 0.0,
                      'max_z_dT_dt': 0.0,
                      'nusselt': [],
                      }
        self.layers = self.layers.append(layers_row, ignore_index=True)
        self.matrix = True
        console.event("Finished inserting matrix ({}) into the model! (task took {}s)".format(material,
                                   time.time() - t_start), verbose=self.verbose)


    def insert_boundary(self, temperature, depth_range, material='boundary'):
        """
        Insert boundary layers into the linear model within a specified depth range.
        :param temperature:
        :param depth_range:
        :param material:
        :return:
        """
        console.event("Inserting boundary ({}) into the box!".format(material), verbose=self.verbose)
        t_start = time.time()
        subdf = self.mesh.query('{} <= coords <= {}'.format(depth_range[0], depth_range[1]))
        coords = subdf['coords'].tolist()
        temperatures = self.mesh['temperature'].tolist()
        objects = self.mesh['object'].tolist()
        object_ids = self.mesh['object_id'].tolist()
        for coord in coords:
            if depth_range[0] <= coord <= depth_range[1]:
                index = backends.predict_linear_index(z=coord, spatial_res=self.spatial_res)
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

    def to_csv(self):
        """
        Write the mesh dataframe to a csv.
        :return:
        """
        console.event("mesh.csv now available in {}!".format(os.getcwd()), verbose=self.verbose)
        self.mesh.to_csv("mesh.csv")

    def get_mesh(self):
        """
        Return the mesh dataframe.
        :return:
        """
        return self.mesh

    def get_layers(self):
        """
        Returns the layers dataframe.
        :return:
        """
        return self.layers

    def verify_box(self):
        """
        Verify that the model has been constructed properly.
        :return:
        """
        objects = self.mesh['object'].tolist()
        if 'NONE' in objects:
            console.error("WARNING! Not all coordinate spaces defined in box!", verbose=True)
        else:
            console.event("Box verified!", verbose=True)

    def update(self, auto_update=True, timestep=False):
        """
        Update the model over one time step. Has the ability to run the model to completion.
        :param auto_update:
        :param timestep:
        :return:
        """
        t = time.time()
        while auto_update is True and self.evolution_time > 0:
            self.delta_time = backends.override_timestep(timestep=timestep, conductivities=self.conductivities,
                                                spatial_res=self.spatial_res, spatial_sigfigs=self.spatial_sigfigs)
            console.nominal("Model time at: {} (timestep: {})...".format(
                self.evolution_time, self.delta_time), verbose=self.verbose)
            z = self.mesh['coords'].tolist()
            z_plus = self.mesh['zplus_index'].tolist()
            z_minus = self.mesh['zminus_index'].tolist()
            object_ids = self.mesh['object_id'].tolist()
            temperatures = self.mesh['temperature'].tolist()
            conductivity = self.mesh['conductivity'].tolist()
            conduction = heat.conduction_linear(z=z, z_plus_indeces=z_plus, z_minus_indeces=z_minus,
                                        temperatures=temperatures, conductivity=conductivity,
                                        spatial_res=self.spatial_res, delta_time=self.delta_time, object_ids=object_ids)
            self.mesh['temperature'] = conduction[0]
            self.mesh['dT_dt'] = conduction[1]
            console.nominal("Finished modeling conduction! (task took {}s)".format(time.time() - t), verbose=self.verbose)
            heat.nusselt(layers=self.layers, mesh=self.mesh, spatial_res=self.spatial_res)
            new_evolution_time = round(self.evolution_time - self.delta_time, self.spatial_sigfigs)
            self.evolution_time = new_evolution_time

        console.event("Model time is at 0! (task took {}s)".format(time.time() - t), verbose=self.verbose)
        return None
