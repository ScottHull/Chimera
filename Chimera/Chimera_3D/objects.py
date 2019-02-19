import numpy as np
from Chimera.Chimera_3D import dynamics, backends, heat, console

def object_actions(objects_df, coords, matrix_ids, spatial_res, spatial_sigfigs, evolution_time, initial_time, max_x, max_y, max_z,
                   x_plus, x_minus, y_plus, y_minus, z_plus, z_minus, delta_time, matrix_densities, matrix_viscosities,
                   mesh_temperatures, matrix_conductivities, lower_model, upper_model, matrix_diffusivities,
                   verbose, chemical_partitioning, chemistry, mesh_fO2, mesh_pressures, mesh_objects, gravity, conduction=True):
    # extract object information & parameters from the objects dataframe contained in the box
    object_objects = np.array(objects_df['object'])
    object_object_ids = np.array(objects_df['object_id'])
    object_temperatures = np.array(objects_df['temperature'])
    object_densities = np.array(objects_df['density'])
    object_coords = np.array(objects_df['coords'], dtype=object)
    object_velocities = np.array(objects_df['velocity'], dtype=object)
    object_conductivities = np.array(objects_df['conductivity'])
    object_radii = np.array(objects_df['radius'])
    cell_vertices = np.array(objects_df['cell_vertices'], dtype=object)
    cell_indices = np.array(objects_df['cell_indices'], dtype=object)
    nearest_indices = np.array(objects_df['nearest_index'], dtype=object)
    drag_coeffs = np.array(objects_df['drag_coeff'])
    compositions = np.array(objects_df['composition'])
    cps = np.array(objects_df['cp'])
    copy_mesh_temperatures = mesh_temperatures
    for object_index, object_object in enumerate(object_objects):
        object_id = object_object_ids[object_index]  # get the current object id
        coord = object_coords[object_index]  # get the current object coordinate
        #  interpolate to find the cell in which the object exists
        #  calculate the object's settling velocity
        cell = backends.interpolate_cell(
            coord=coord,
            spatial_res=spatial_res,
            spatial_sigfigs=spatial_sigfigs,
            max_x=max_x,
            max_y=max_y,
            max_z=max_z,
            x_plus=x_plus,
            x_minus=x_minus,
            y_plus=y_plus,
            y_minus=y_minus,
            z_plus=z_plus,
            z_minus=z_minus,
            verbose=verbose
        )
        nearest_indices[object_index] = cell[0]
        cell_indices[object_index] = cell[1]
        cell_vertices[object_index] = cell[2]
        velocity = dynamics.stokes_terminal(
            density=object_densities[object_index],
            density_matrix=matrix_densities[nearest_indices[object_index]],
            drag_coeff=drag_coeffs[object_index],
            radius=object_radii[object_index],
            viscosity_matrix=matrix_viscosities[nearest_indices[object_index]],
            current_coord=coord,
            lower_model=lower_model,
            upper_model=upper_model,
            gravity=gravity,
        )
        #  get the object's new coordinates based on the object's velocity
        updated_coords, distance_travelled = backends.update_position(
            coord=coord,
            velocity=velocity,
            delta_time=delta_time,
            max_z=lower_model,
            spatial_res=spatial_res
        )
        cell = backends.interpolate_cell(
            coord=updated_coords,
            spatial_res=spatial_res,
            spatial_sigfigs=spatial_sigfigs,
            max_x=max_x,
            max_y=max_y,
            max_z=max_z,
            x_plus=x_plus,
            x_minus=x_minus,
            y_plus=y_plus,
            y_minus=y_minus,
            z_plus=z_plus,
            z_minus=z_minus,
            verbose=verbose
        )
        # determine the mesh cell within which the object resides
        # determine how much heat was produced through drag on the object

        if conduction:
            viscous_heat = heat.viscous_dissipation(
                drag_coeff=drag_coeffs[object_index],
                cp=cps[object_index],
                delta_time=delta_time,
                matrix_density=matrix_densities[cell[0]],
                object_density=object_densities[object_index],
                object_radius=object_radii[object_index],
                object_velocity=velocity,
                distance_travelled=distance_travelled
            )
            # add the heat generated from viscous dissipation to the object
            object_temperatures[object_index] += viscous_heat

        # allow the heat from the object to conduct, assuming only in contact with nearest cell vertex
        if conduction:
            heat.object_conduction(
                object_temperatures=object_temperatures,
                copy_mesh_temperatures=copy_mesh_temperatures,
                object_index=object_index,
                object_k=object_conductivities[object_index],
                spatial_res=spatial_res,
                delta_time=delta_time,
                mesh_temperatures=mesh_temperatures,
                nearest_index=cell[0],
                farthest_index=cell[3],
                directional_vertices=cell[4],
                vertex_distances=cell[5],
                total_distance=cell[6],
                matrix_k=matrix_conductivities[cell[0]],
                matrix_ids=matrix_ids,
                coords=coords,
                matrix_diffusivities=matrix_diffusivities,
                spatial_sigfigs=spatial_sigfigs,
                verbose=verbose
            )
        # chemically equilibrate the object with its surroundings
        if chemical_partitioning:
            chemistry.equilibrate(
                object_moles=compositions,
                object_index=object_index,
                vertex_distances=cell[5],
                matrix_ids=matrix_ids,
                total_distance=cell[6],
                vertex_indices=cell[1],
                pressures=mesh_pressures,
                temperatures=mesh_temperatures,
                object_temperatures=object_temperatures,
                fO2=mesh_fO2,
                spatial_res=spatial_res,
                object_radius=object_radii[object_index],
                object_id=object_id,
                z_depth=coord[2],
            )
        # console.event("{} ({}) will travel from {} to {} (velocity: {})".format(
        #     object_object, object_id, coord, updated_coords, velocity), verbose=verbose)

        #  update the dataframes with the new data
        object_velocities[object_index] = velocity
        object_coords[object_index] = updated_coords
        nearest_indices[object_index] = cell[0]  # nearest index
        cell_indices[object_index] = cell[1]  # cell indices
        cell_vertices[object_index] = cell[2]  # cell vertices
    objects_df['coords'] = object_coords
    objects_df['velocity'] = object_velocities
    objects_df['cell_vertices'] = cell_vertices
    objects_df['cell_indices'] = cell_indices
    objects_df['nearest_index'] = nearest_indices
    objects_df['temperature'] = object_temperatures
    return object_coords, nearest_indices, cell_indices