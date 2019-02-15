import sys
from functools import partial
import pandas as pd
import numpy as np
import time
import multiprocessing as mp
from random import uniform
from Chimera.Chimera_3D import mesh, console, backends, neighbors, heat, plots, objects, chemistry, mainloop
import warnings; warnings.filterwarnings('ignore')
# import pyximport; pyximport.install()


class Box:

    def __init__(self, evolution_time, gravity, multiprocessing=False, num_processors=1, conduction=True,
                 settling_mode='stokes terminal', radioactivity=True, chem=True, verbose=True):
        self.mesh = pd.DataFrame(
            {
        }
        )
        self.objects = pd.DataFrame(
            {
        }
        )
        self.gravity = gravity,
        self.conduction = conduction
        self.settling_mode = settling_mode
        self.radioactivity = radioactivity
        self.initial_time = evolution_time
        self.evolution_time = evolution_time
        self.verbose = verbose
        self.max_x = 0.0
        self.max_y = 0.0
        self.max_z = 0.0
        self.spatial_res = 0
        self.spatial_sigfigs = 0.0
        self.dimension = None
        self.matrix = None
        self.object = None
        self.boundary = None
        self.chem = chem
        self.chemicaldiffusion = None
        self.lower_model = self.max_z
        self.upper_model = 0.0
        self.id_val = 0
        self.conductivities = []
        self.plots = plots.plots()
        self.iterations = 0
        self.multiprocessing = multiprocessing
        self.num_workers = num_processors
        if self.num_workers > mp.cpu_count():
            self.num_workers = mp.cpu_count()
        self.chemistry = chemistry.Chemistry(box=self)
        self.objectGenerators = {

        }
        self.num_objects = 0
        self.log_interval = None
        self.num_nodes = 0


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
        self.spatial_sigfigs = int(m.get_spatial_sigfigs())
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
        neighbor_count = 1
        total_count = len(self.mesh['coords'])
        t_start = time.time()
        arr = np.array(self.mesh['coords'])  # an array of all the mesh coordinates
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

    def insert_matrix(self, material, temperature, conductivity, fO2, density, viscosity, depth_range, heat_capacity,
                      pressure, composition={}, element_diffusivity={}, grad_temperature=0.0,
                      grad_pressure=0.0, grad_fO2=0.0):
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
        coords = np.array(subdf['coords'])
        temperatures = np.array(self.mesh['temperature'])
        objects = np.array(self.mesh['object'])
        object_ids = np.array(self.mesh['object_id'])
        conductivities = np.array(self.mesh['conductivity'])
        densities = np.array(self.mesh['density'])
        viscosities = np.array(self.mesh['viscosity'])
        diffusivities = np.array(self.mesh['diffusivity'])
        fO2s = np.array(self.mesh['fO2'])
        pressures = np.array(self.mesh['pressure'])
        self.conductivities.append(conductivity)
        len_coords = len(coords)
        # insert the matrix into coordinate positions defined by the user's z-range
        # TODO: we can limit recursion time by using the index prediction equation here
        for coord in coords:
            # z = coords[2][self.dimension - 1]
            z = coord[2]
            if depth_range[0] <= z <= depth_range[1]:
                index = backends.predict_index(
                    coord=coord,
                    max_x=self.max_x,
                    max_y=self.max_y,
                    max_z=self.max_z,
                    spatial_res=self.spatial_res,
                    verbose=self.verbose
                )
                console.nominal("Inserting matrix ({}) at {}...".format(material, coord), verbose=self.verbose)
                conductivities[index] = float(conductivity)
                diffusivities[index] = float(conductivity / (density * heat_capacity))
                objects[index] = material
                densities[index] = float(density)
                viscosities[index] = float(viscosity)
                temperatures[index] = float(temperature + (grad_temperature * round(((z -
                                    depth_range[0]) / self.spatial_res), self.spatial_sigfigs)))
                pressures[index] = float(pressure + (grad_pressure * round(((z -
                                    depth_range[0]) / self.spatial_res), self.spatial_sigfigs)))
                fO2s[index] = float(fO2 + (grad_fO2 * round(((z -
                                    depth_range[0]) / self.spatial_res), self.spatial_sigfigs)))
                object_ids[index] = backends.generate_object_id(object_type='matrix',
                                                                id_val=self.id_val)
                if self.chem:
                    self.chemistry.insertMatrixComposition(index=index, composition=composition, material=material,
                                                           diffusivity=element_diffusivity)
                self.id_val += 1
        # reset the dataframe columns
        self.mesh['temperature'] = temperatures
        self.mesh['object'] = objects
        self.mesh['object_id'] = object_ids
        self.mesh['density'] = densities
        self.mesh['viscosity'] = viscosities
        self.mesh['conductivity'] = conductivities
        self.mesh['diffusivity'] = diffusivities
        self.mesh['fO2'] = fO2s
        self.mesh['pressure'] = pressures
        self.matrix = True
        # insert composition into the chemistry instance so it can be tracked independently to avoid
        # box getting garbled up
        console.event("Finished inserting matrix ({}) into the box! (task took {}s)".format(material,
                                time.time() - t_start), verbose=self.verbose)
        return material

    def insert_boundary(self, temperature, depth_range, location, material='boundary', composition={},
                        element_diffusivity={}):
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
        coords = np.array(subdf['coords'])
        temperatures = np.array(self.mesh['temperature'])
        objects = np.array(self.mesh['object'])
        object_ids = np.array(self.mesh['object_id'])
        for coord in coords:
            z = coord[2]
            if depth_range[0] <= z <= depth_range[1]:
                index = backends.predict_index(
                    coord=coord,
                    max_x=self.max_x,
                    max_y=self.max_y,
                    max_z=self.max_z,
                    spatial_res=self.spatial_res,
                    verbose=self.verbose
                )
                console.nominal("Inserting boundary ({}) at {}...".format(material, coord), verbose=self.verbose)
                temperatures[index] = temperature
                objects[index] = material
                object_ids[index] = backends.generate_object_id(
                    object_type='boundary',
                    id_val=self.id_val
                )
                if self.chem:
                    self.chemistry.insertMatrixComposition(index=index, composition=composition, material=material,
                                                           diffusivity=element_diffusivity)
                self.id_val += 1
        self.mesh['temperature'] = temperatures
        self.mesh['object'] = objects
        self.mesh['object_id'] = object_ids
        self.boundary = True
        # the upper-most usable portion of the model is the matrix immediately below the boundary
        if location.lower() == "top" or location.lower() == "t":
            self.upper_model = 0.0 + self.spatial_res
        elif location.lower() == "bottom" or location.lower() == "b":
            self.lower_model = self.max_z - self.spatial_res
        console.event("Finished inserting boundary ({}) into the box! (task took {}s)".format(
            material, time.time() - t_start), verbose=self.verbose
        )

        return

    def insert_object(self, material, temperature, radius, conductivity, density, drag_coeff, cp, x, y, z,
                      composition={}):
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
        drag_coeffs = self.objects['drag_coeff'].tolist()
        cps = self.objects['cp'].tolist()
        compositions = self.objects['composition'].tolist()
        console.nominal(
            "Inserting object ({}) at {}...".format(
                material, requested_coord), verbose=self.verbose
        )
        objects.append(material)
        object_ids.append(
            backends.generate_object_id(
                object_type='object',
                id_val=self.id_val
            )
        )
        radii.append(float(radius))
        locs.append(requested_coord)
        temperatures.append(float(temperature))
        conductivities.append(float(conductivity))
        densities.append(float(density))
        velocities.append((0.0, 0.0, 0.0))
        drag_coeffs.append(float(drag_coeff))
        cps.append(float(cp))
        compositions.append(composition)
        # object_headers = [self.objects.columns.values]
        self.objects = pd.DataFrame()  # reset the dataframe to avoid length issues
        self.objects['object'] = objects
        self.objects['object_id'] = object_ids
        self.objects['temperature'] = temperatures
        self.objects['radius'] = radii
        self.objects['coords'] = locs
        self.objects['conductivity'] = conductivities
        self.objects['density'] = densities
        self.objects['drag_coeff'] = drag_coeffs
        self.objects['cp'] = cps
        self.objects['velocity'] = velocities
        self.objects['nearest_index'] = np.NAN
        self.objects['cell_indices'] = np.NAN
        self.objects['cell_vertices'] = np.NAN
        self.objects['composition'] = compositions
        self.id_val += 1
        self.object = True
        self.num_objects += 1
        console.event(
            "Finished inserting object ({}) into the box! (task took {}s)".format(
                material, time.time() - t_start), verbose=self.verbose
        )

        return

    def insertObjectGenerator(self, generator_name, material, num_per_iteration, temperature, radius, conductivity,
                              density, drag_coeff, cp, x_range, y_range, z_range, composition={}):

        def objectGenerator(instance, material, num_per_iteration, temperature, radius, conductivity,
                            density, drag_coeff, cp, x_range, y_range, z_range, composition={}):

            for i in range(num_per_iteration):
                x = uniform(x_range[0], x_range[1])
                y = uniform(y_range[0], y_range[1])
                z = uniform(z_range[0], z_range[1])
                instance.insert_object(
                    material=material,
                    conductivity=conductivity,
                    temperature=temperature,
                    radius=radius,
                    density=density,
                    drag_coeff=drag_coeff,
                    cp=cp,
                    x=x,
                    y=y,
                    z=z,
                    composition=composition
                )

        self.objectGenerators.update(
            {
                generator_name: partial(
                    objectGenerator,
                    instance=self,
                    material=material,
                    num_per_iteration=num_per_iteration,
                    temperature=temperature,
                    radius=radius,
                    conductivity=conductivity,
                    density=density,
                    drag_coeff=drag_coeff,
                    cp=cp,
                    x_range=x_range,
                    y_range=y_range,
                    z_range=z_range,
                    composition=composition
                )
            }
        )




    def verify_box(self):
        """
        Verifies that the box is constructed properly.
        :return:
        """
        objects = np.array(self.mesh['object'])
        if 'NONE' in objects:
            console.error("WARNING! Not all coordinate spaces defined in box!", verbose=True)
        else:
            console.event("Box verified!", verbose=self.verbose)

    def to_csv(self, model_time=""):
        """
        Returns csv formats of all stored dataframes.
        :return:
        """
        self.mesh.to_csv("mesh_{}.csv".format(model_time), index=False)
        self.objects.to_csv("objects_{}.csv".format(model_time), index=False)
        pd.DataFrame(self.chemistry.track_distribution_coeffs).to_csv("D_coeffs.csv")
        # matrix_chem_file = open("matrix_chem_{}.txt".format(model_time), 'w')
        # matrix_chem_file.write(str(self.chemistry.matrix))
        # matrix_chem_file.close()
        partition_file = open("partitioning_{}.txt".format(model_time), 'w')
        partition_file.write(str(self.chemistry.partitioning))
        partition_file.close()
        # pd.DataFrame(self.chemistry.matrix).to_csv("matrix_chem_{}.csv".format(model_time), index=False)
        # pd.DataFrame(self.chemistry.partitioning).to_csv("partitioning_{}.csv".format(model_time), index=False)

    def return_mesh(self, mesh, **kwargs):
        for key, val in kwargs.items():
            mesh[str(key)] = val
        return mesh


    def update(self, auto_update=True, timestep=False, animate_model=False, show_model=False, log_interval=None):
        """
        Update the model over one time step. Has the ability to run the model to completion.
        :param auto_update:
        :param timestep:
        :return:
        """
        t = time.time()
        # extract values that are not updated during the loop from the dataframe as lists
        coords = np.array(self.mesh['coords'])
        x_plus = np.array(self.mesh['xplus_index'])
        x_minus = np.array(self.mesh['xminus_index'])
        y_plus = np.array(self.mesh['yplus_index'])
        y_minus = np.array(self.mesh['yminus_index'])
        z_plus = np.array(self.mesh['zplus_index'])
        z_minus = np.array(self.mesh['zminus_index'])
        object_ids = np.array(self.mesh['object_id'])
        mesh_objects = np.array(self.mesh['object'])
        conductivities = np.array(self.mesh['conductivity'])
        viscosity = np.array(self.mesh['viscosity'])
        density = np.array(self.mesh['density'])
        diffusivities = np.array(self.mesh['diffusivity'])
        temperatures = np.array(self.mesh['temperature'])  # load in current temperatures across the mesh
        fO2s = np.array(self.mesh['fO2'])
        pressures = np.array(self.mesh['pressure'])
        dT_dts = np.array(self.mesh['dT_dt'])
        mesh_indices = np.array(self.mesh.index)
        len_coords = len(coords)
        # calculate the timestep based on the maximum conductivity of material in the box
        self.delta_time = backends.override_timestep(
            timestep=timestep,
            conductivities=conductivities,
            spatial_res=self.spatial_res,
            spatial_sigfigs=self.spatial_sigfigs,
            diffusivities=diffusivities,
            verbose=self.verbose
        )
        # will run model to completion while the remaining time is above 0
        while auto_update is True and self.evolution_time > 0:
            console.nominal("Model time at: {} (timestep: {})...".format(
                self.evolution_time, self.delta_time), verbose=self.verbose)
            # run the object generator functions stored in the instance, if any
            for i in self.objectGenerators:
                self.objectGenerators[i]()
            # temperatures = np.array(self.mesh['temperature'])
            #  load in current temperatures across the mesh
            #  perform actions on objects inside of the model but independent of the mesh
            object_coords, nearest_indices, cell_indices = objects.object_actions(
                mesh_temperatures=temperatures,
                objects_df=self.objects,
                spatial_res=self.spatial_res,
                spatial_sigfigs=self.spatial_sigfigs,
                evolution_time=self.evolution_time,
                delta_time=self.delta_time,
                initial_time=self.initial_time,
                matrix_densities=density,
                matrix_viscosities=viscosity,
                x_plus=x_plus,
                x_minus=x_minus,
                y_plus=y_plus,
                y_minus=y_minus,
                z_plus=z_plus,
                z_minus=z_minus,
                max_x=self.max_x,
                max_y=self.max_y,
                max_z=self.max_z,
                lower_model=self.lower_model,
                upper_model=self.upper_model,
                verbose=self.verbose,
                conduction=self.conduction,
                matrix_conductivities=conductivities,
                matrix_ids=object_ids,
                coords=coords,
                matrix_diffusivities=diffusivities,
                chem=self.chem,
                chemistry=self.chemistry,
                mesh_fO2=fO2s,
                mesh_objects=mesh_objects,
                mesh_pressures=pressures,
                gravity=self.gravity,
            )
            # if a PDE has to be calculated, loop through the model
            # will pass if no PDE is needed in order to optimize runtime
            if self.conduction or self.chem:
                loop_t = time.time()
                loop_outputs = mainloop.modelLoop(
                    conduction=self.conduction,
                    chem=self.chem,
                    coords=coords,
                    chemistry=self.chemistry,
                    len_coords=len_coords,
                    x_plus_indices=x_plus,
                    x_minus_indices=x_minus,
                    y_plus_indices=y_plus,
                    y_minus_indices=y_minus,
                    z_plus_indices=z_plus,
                    z_minus_indices=z_minus,
                    temperatures=temperatures,
                    object_ids=object_ids,
                    spatial_res=self.spatial_res,
                    conductivities=conductivities,
                    delta_time=self.delta_time,
                    mesh_indices=mesh_indices,
                    num_workers=self.num_workers,
                    spatial_sigfigs=self.spatial_sigfigs,
                    therm_diffusivities=diffusivities,
                    verbose=self.verbose,
                    multiprocess=False
                )
                if self.conduction:
                    temperatures = loop_outputs[0]
                    dT_dts = loop_outputs[1]
                if self.chem:
                    # self.chemistry.matrix = loop_outputs[2]
                    self.chemistry.resetMatrixComp(new_matrix_comp=loop_outputs[2])
            # plot the model's dynamic components, remove in the operational version
            self.plots.plot_cell_therm(
                object_coords=object_coords,
                nearest_coords=nearest_indices,
                vertex_indices=cell_indices,
                mesh_coords=coords,
                max_x=self.max_x,
                max_y=self.max_y,
                max_z=self.max_z,
                temperatures=temperatures,
                spatial_res=self.spatial_res,
                model_time=self.evolution_time,
                save=animate_model,
                show=show_model,
                heat=self.conduction
            )
            self.plots.plot_cell_chem(
                object_coords=object_coords,
                vertex_indices=cell_indices,
                mesh_coords=coords,
                max_x=self.max_x,
                max_y=self.max_y,
                max_z=self.max_z,
                mesh_composition=self.chemistry.matrix,
                element='w',
                model_time=self.evolution_time,
                spatial_res=self.spatial_res,
                chem=self.chem,
                save=animate_model,
                show=show_model
            )
            # create logs at model intervals
            if log_interval is not None:
                if self.log_interval == log_interval:
                    self.to_csv(model_time=self.evolution_time)
                    self.log_interval = 0
                else:
                    self.log_interval += 1
            # update the new time in the model
            new_evolution_time = round(self.evolution_time - self.delta_time, self.spatial_sigfigs)
            self.evolution_time = new_evolution_time
            self.iterations += 1

        console.event(
            "Model time is at 0! (task took {}s ({} iterations (1 iteration = {}s), {} timestep)".format(
            time.time() - t, self.iterations, (time.time() - t) / self.iterations, self.delta_time),
            verbose=self.verbose
        )
        if self.chem:
            self.mesh['composition'] = np.array(self.chemistry.matrix)
            self.return_mesh(
                mesh=self.mesh,
                temperature=temperatures,
                dT_dt=dT_dts,
                conductivity=conductivities,
                density=density,
                viscosity=viscosity,
                mesh_composition=np.array(self.chemistry.matrix)
            )
        else:
            self.return_mesh(
                mesh=self.mesh,
                temperature=temperatures,
                dT_dt=dT_dts,
                conductivity=conductivities,
                density=density,
                viscosity=viscosity,
            )
        # will create animations of models if specified
        if animate_model is True:
            self.plots.animate(initial_time=self.initial_time, chem=self.chem, conduction=self.conduction)
        return None
