import sys
import pandas as pd
import numpy as np
import time
from Chimera.Chimera_3D import mesh, console, backends, neighbors, heat, plots, settling_modes
import warnings; warnings.filterwarnings('ignore')
# import pyximport; pyximport.install()


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
        self.conductivities = []
        self.plots = plots.plots()

    def build(self, spatial_res, x, y, z=None):
        """
        Constructs the mesh based on specified parameters.
        :param spatial_res:
        :param x:
        :param y:
        :param z:
        :return:
        """
        if z is not 3:
            self.dimension = 3
            self.max_x, self.max_y, self.max_z = x, y, z
        else:
            self.dimension = 2
            self.max_x, self.max_y, self.max_z = x, y, None
        self.spatial_res = spatial_res
        console.event("Constructing box...", verbose=self.verbose)
        console.nominal("Building mesh...", verbose=self.verbose)
        # instantiate the mesh
        m = mesh.Mesh(spatial_res=spatial_res, x=x, y=y, z=z, verbose=self.verbose)
        self.spatial_sigfigs = m.get_spatial_sigfigs()
        console.nominal("Assigning mesh to dataframe...", verbose=self.verbose)
        # construct the mesh
        m.build(df=self.mesh)
        console.nominal("Exiting box construction...", verbose=self.verbose)
        console.event("Box constructed!", verbose=self.verbose)
        console.event("Fetching nearest neighbors...", verbose=self.verbose)
        coords_range = range(len(self.mesh['coords']))
        # create all of the columns for the instantiated dataframes
        backends.set_columns(mesh_df=self.mesh, object_df=self.objects, coords_range=coords_range)
        self.x_coords = np.NAN
        self.y_coords = np.NAN
        self.z_coords = np.NAN
        spatial_sigfigs = m.get_spatial_sigfigs()  # find the number of sigfigs defined by the spatial resolution
        neighbor_count = 0
        total_count = len(self.mesh['coords'])
        t_start = time.time()
        arr = self.mesh['coords'].tolist()  # an array of all the mesh coordinates
        x_plus = []
        x_minus = []
        y_plus = []
        y_minus = []
        z_plus = []
        z_minus = []
        # find nearest neighbors in 3D
        for coords in arr:
            if self.dimension is 3:
                console.nominal("Finding neighbor for ({}, {}, {}) ({}/{} points)".format(
                    coords[0], coords[1], coords[2], neighbor_count, total_count), verbose=self.verbose)
            else:
                console.nominal("Finding neighbor for ({}, {})".format(
                    coords[0], coords[1]), verbose=self.verbose)
            n = neighbors.get_neighbors(verbose=self.verbose, coords=coords, max_x=x, max_y=y, max_z=z,
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

    def insert_matrix(self, material, temperature, conductivity, density, viscosity, depth_range):
        """
        Insert a matrix into the model.
        :param material:
        :param temperature:
        :param conductivity:
        :param depth_range:
        :return:
        """
        console.event("Inserting matrix ({}) into the box!".format(material), verbose=self.verbose)
        t_start = time.time()
        subdf = self.mesh.query('{} <= z_vals <= {}'.format(depth_range[0], depth_range[1]))
        self.mesh.drop(['z_vals'], axis=1)
        coords = subdf['coords'].tolist()
        temperatures = self.mesh['temperature'].tolist()
        objects = self.mesh['object'].tolist()
        object_ids = self.mesh['object_id'].tolist()
        conductivities = self.mesh['conductivity'].tolist()
        densities = self.mesh['density'].tolist()
        viscosities = self.mesh['viscosity'].tolist()
        self.conductivities.append(conductivity)
        len_coords = len(coords)
        # insert the matrix into coordinate positions defined by the user's z-range
        # TODO: we can limit recursion time by using the index prediction equation here
        for coord in coords:
            if depth_range[0] <= coords[2][self.dimension - 1] <= depth_range[1]:
                index = backends.predict_index(coord=coord, max_x=self.max_x, max_y=self.max_y, max_z=self.max_z,
                                               spatial_res=self.spatial_res, verbose=self.verbose)
                console.nominal("Inserting matrix ({}) at {}...".format(material, coord), verbose=self.verbose)
                temperatures[index] = temperature
                conductivities[index] = conductivity
                objects[index] = material
                densities[index] = density
                viscosities[index] = viscosity
                object_ids[index] = backends.generate_object_id(object_type='matrix',
                                                                id_val=self.id_val)
                self.id_val += 1
        self.mesh['temperature'] = temperatures
        self.mesh['object'] = objects
        self.mesh['object_id'] = object_ids
        self.mesh['density'] = densities
        self.mesh['viscosity'] = viscosities
        self.mesh['conductivity'] = conductivities
        self.matrix = True
        console.event("Finished inserting matrix ({}) into the box! (task took {}s)".format(material,
                                time.time() - t_start), verbose=self.verbose)

    def insert_boundary(self, temperature, depth_range, material='boundary'):
        """
        Insert boundary layers in the model.
        :param temperature:
        :param depth_range:
        :param material:
        :return:
        """
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

    def insert_object(self, material, temperature, radius, conductivity, density, x, y, z=None):
        """
        Insert an object in the the objects dataframe.
        :param material:
        :param temperature:
        :param radius:
        :param conductivity:
        :param x:
        :param y:
        :param z:
        :return:
        """
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
        locs = self.objects['coords'].tolist()
        densities = self.objects['density'].tolist()
        velocities = self.objects['velocity'].tolist()
        conductivities = self.objects['conductivity'].tolist()
        console.nominal("Inserting object ({}) at {}...".format(material, requested_coord), verbose=self.verbose)
        objects.append(material)
        object_ids.append(backends.generate_object_id(object_type='object',
                                                      id_val=self.id_val))
        radii.append(radius)
        locs.append(requested_coord)
        temperatures.append(temperature)
        conductivities.append(conductivity)
        densities.append(density)
        velocities.append(0.0)
        object_headers = [self.objects.columns.values]
        self.objects = pd.DataFrame()  # reset the dataframe to avoid length issues
        self.objects['object'] = objects
        self.objects['object_id'] = object_ids
        self.objects['temperature'] = temperatures
        self.objects['radius'] = radii
        self.objects['coords'] = locs
        self.objects['conductivity'] = conductivities
        self.objects['density'] = densities
        self.objects['velocity'] = velocities
        self.objects['nearest_index'] = np.NAN
        self.objects['cell_indices'] = np.NAN
        self.objects['cell_vertices'] = np.NAN
        self.id_val += 1
        self.object = True
        console.event("Finished inserting object ({}) into the box! (task took {}s)".format(material,
                                                                        time.time() - t_start), verbose=self.verbose)

    def verify_box(self):
        """
        Verifies that the box is constructed properly.
        :return:
        """
        objects = self.mesh['object'].tolist()
        if 'NONE' in objects:
            console.error("WARNING! Not all coordinate spaces defined in box!", verbose=True)
            self.mesh.drop(['z_vals'])
        else:
            console.event("Box verified!", verbose=True)

    def to_csv(self):
        """
        Returns csv formats of all stored dataframes.
        :return:
        """
        self.mesh.to_csv("mesh.csv", index=False)
        self.objects.to_csv("objects.csv", index=False)

    def update(self, auto_update=True, timestep=False, animate_model=False):
        """
        Update the model over one time step. Has the ability to run the model to completion.
        :param auto_update:
        :param timestep:
        :return:
        """
        t = time.time()
        # extract values that are not updated during the loop from the dataframe as lists
        coords = self.mesh['coords'].tolist()
        x_plus = self.mesh['xplus_index'].tolist()
        x_minus = self.mesh['xminus_index'].tolist()
        y_plus = self.mesh['yplus_index']
        y_minus = self.mesh['yminus_index']
        z_plus = self.mesh['zplus_index'].tolist()
        z_minus = self.mesh['zminus_index'].tolist()
        object_ids = self.mesh['object_id'].tolist()
        conductivity = self.mesh['conductivity'].tolist()
        viscosity = self.mesh['viscosity'].tolist()
        density = self.mesh['density'].tolist()
        # calculate the timestep based on the maximum conductivity of material in the box
        self.delta_time = backends.override_timestep(timestep=timestep, conductivities=self.conductivities,
                                                     spatial_res=self.spatial_res, spatial_sigfigs=self.spatial_sigfigs)
        # will run model to completion while the remaining time is above 0
        while auto_update is True and self.evolution_time > 0:
            console.nominal("Model time at: {} (timestep: {})...".format(
                self.evolution_time, self.delta_time), verbose=self.verbose)
            temperatures = self.mesh['temperature'].tolist()
            #  perform actions on objects inside of the model but independent of the mesh
            object_coords, nearest_indices, cell_indices = backends.object_actions(mesh_df=self.mesh, objects_df=self.objects,
                                    spatial_res=self.spatial_res, spatial_sigfigs=self.spatial_sigfigs,
                                    evolution_time=self.evolution_time, delta_time=self.delta_time,
                                    initial_time=self.initial_time, matrix_densities=density,
                                    matrix_viscosities=viscosity, x_plus=x_plus, x_minus=x_minus, y_plus=y_plus,
                                    y_minus=y_minus, z_plus=z_plus, z_minus=z_minus, max_x=self.max_x, max_y=self.max_y,
                                    max_z=self.max_z, verbose=self.verbose)
            # plot the model's dynamic components
            self.plots.plot_cell(object_coords=object_coords, nearest_coords=nearest_indices,
                            vertex_indeces=cell_indices, mesh_coords=coords, max_x=self.max_x, max_y=self.max_y,
                            max_z=self.max_z, spatial_res=self.spatial_res, model_time=self.evolution_time,
                            save=animate_model, show=False)
            # finite central difference conductivity across entire box
            conduction = heat.conduction(coords=coords, x_plus_indeces=x_plus, x_minus_indeces=x_minus,
                                        y_plus_indeces=y_plus, y_minus_indeces=y_minus, z_plus_indeces=z_plus,
                                        z_minus_indeces=z_minus, temperatures=temperatures, conductivity=conductivity,
                                        spatial_res=self.spatial_res, delta_time=self.delta_time, object_ids=object_ids)
            # conduction will return tuple: temperature at index 0 and dT/dt at index 1
            self.mesh['temperature'] = conduction[0]
            self.mesh['dT_dt'] = conduction[1]
            console.nominal("Finished modeling conduction! (task took {}s)".format(time.time() - t), verbose=self.verbose)
            # update the new time in the model
            new_evolution_time = round(self.evolution_time - self.delta_time, self.spatial_sigfigs)
            self.evolution_time = new_evolution_time

        console.event("Model time is at 0! (task took {}s)".format(time.time() - t), verbose=self.verbose)
        if animate_model is True:
            self.plots.animate(initial_time=self.initial_time)
        return None
