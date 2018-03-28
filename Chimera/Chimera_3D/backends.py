from Chimera.Chimera_3D import console
import numpy as np
from . import settling_modes

def set_columns(mesh_df, object_df, coords_range):

    mesh_columns = ["object", "object_id", "xplus_index", "xminus_index", "yplus_index", "yminus_index",
                    "zplus_index", "zminus_index", "temperature", "conductivity", "viscosity", "density"]
    object_columns = ["object", "object_id", "coords", "temperature", "radius", "conductivity", "density", "velocity",
                      "cell_indices", "cell_vertices", "nearest_index"]
    for i in mesh_columns:
        if i == "object" or i == "object_id":
            mesh_df[i] = ["NONE" for i in coords_range]
        elif "index" in str(i):
            mesh_df[i] = [0 for i in coords_range]
        elif i == "coords" or i == "velocity":
            mesh_df[i] = [(0,0,0) for i in coords_range]
        else:
            mesh_df[i] = [0.0 for i in coords_range]
    for i in object_columns:
        object_df[i] = np.NAN
    return mesh_columns, mesh_df

def generate_object_id(object_type, id_val):
    """
    Generates object ID codes so that specific objects and materials can be tracked
    Object/matrial types are unique and coded by the a letter specifying the general type followed by a
    unique number combination
    The general object/material types are coded by the first letter as follows:
        - 'A' = object
        - 'B' = matrix
        - 'C' = boundary
        - 'Z' = else
    :param matrix:
    :return: object_id
    """

    def random_gen(object_identifier, id_val):
        object_id = object_identifier + str(id_val)
        return object_id

    object_type = object_type.lower()

    if object_type is 'matrix':
        object_id = random_gen(object_identifier='B', id_val=id_val)
        # while object_id in object_ids:
        #     object_id = random_gen(object_identifier='B')
        return object_id
    elif object_type is 'object':
        object_id = random_gen(object_identifier='A', id_val=id_val)
        return object_id
    elif object_type is 'boundary':
        object_id = random_gen(object_identifier='C', id_val=id_val)
        return object_id
    else:
        object_id = random_gen(object_identifier='Z', id_val=id_val)
        return object_id

def predict_index(coord, max_x, max_y, spatial_res, max_z=None, verbose=True):
    """
    Predicts the index position of a coordinate in the mesh within the mesh dataframe.
    :param coord:
    :param max_x:
    :param max_y:
    :param spatial_res:
    :param max_z:
    :param verbose:
    :return:
    """
    i, j, k = coord[0], coord[1], coord[2]
    s = spatial_res
    x, y, z = max_x, max_y, max_z
    if z is not None:
        index = int(round((((i*y*z)/(s**3)) + (((i*y) + (i*z) + (j*z))/(s**2)) + ((i+j)/s) + (k/s))))
        return index
    else:
        console.error("2D point prediction not implemented!", verbose=verbose)

def override_timestep(timestep, conductivities, spatial_res, spatial_sigfigs):
    """
    Sets the time intervals for model evolution.  An override allows the user to define a custom time interval
    with the risk of poor convergence.
    :param override:
    :return:
    """
    if timestep is False:
        # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        # the explicit method timestep is only stable if:
        # deltaT <= (deltaX^2)/(2*(conductivity))
        # stability requires the maximum usage of conductivity to get the appropriate timestep
        # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        max_conductivity = max(conductivities)
        delta_time = (spatial_res**2) / (2 * (max_conductivity))
        delta_time = round(delta_time, spatial_sigfigs)
        return delta_time
    else:
        delta_time = timestep
        return delta_time

def interpolate_cell(coord, spatial_sigfigs, spatial_res, max_x, max_y, max_z, x_plus, x_minus, y_plus, y_minus,
                     z_plus, z_minus, verbose=True):
    """
    Takes a given coordinate and determines within which cell in the mesh it resides.
    :param coord:
    :param spatial_sigfigs:
    :param spatial_res:
    :param max_x:
    :param max_y:
    :param max_z:
    :param x_plus:
    :param x_minus:
    :param y_plus:
    :param y_minus:
    :param z_plus:
    :param z_minus:
    :param verbose:
    :return:
    """
    def arbitrary_round(coord, spatial_res, spatial_sigfigs):
        """
        Rounds the given coordinate to coordinates defined in the mesh.
        :param coord:
        :param spatial_res:
        :param spatial_sigfigs:
        :return:
        """
        x, y, z = coord[0], coord[1], coord[2]
        nearest_x = round(round(x / spatial_res) * spatial_res, spatial_sigfigs)
        nearest_y = round(round(y / spatial_res) * spatial_res, spatial_sigfigs)
        nearest_z = round(round(z / spatial_res) * spatial_res, spatial_sigfigs)
        return nearest_x, nearest_y, nearest_z

    def get_neighbors(distance, nearest_coord, spatial_res, spatial_sigfigs):
        """
        Determines if the nearest neighbor should be a "plus coord" or a "minus coord" in space in the given dimensional component.
        :param distance:
        :param nearest_coord:
        :param spatial_res:
        :param spatial_sigfigs:
        :return:
        """
        if distance < 0:
            minus_coord = round(nearest_coord - spatial_res, spatial_sigfigs)
            return minus_coord
        elif distance > 0:
            plus_coord = round(nearest_coord + spatial_res, spatial_sigfigs)
            return plus_coord
        else:
            return nearest_coord

    x, y, z = coord[0], coord[1], coord[2]  # the x, y, z components of the coordinate
    nearest_x, nearest_y, nearest_z = arbitrary_round(coord=coord, spatial_res=spatial_res,
                                spatial_sigfigs=spatial_sigfigs)  # find the nearest defined coordinates in the mesh
    nearest_coord = (nearest_x, nearest_y, nearest_z)  # the nearest coordinate point (tuple)
    nearest_index = predict_index(coord=nearest_coord, max_x=max_x, max_y=max_y,
                    max_z=max_z, spatial_res=spatial_res, verbose=verbose)  # get the mesh index position of the nearest coordinate point
    distance_x, distance_y, distance_z = (x - nearest_x), (y - nearest_y), (z - nearest_z)  # find the distance of the real coordinate to the nearest mesh coordinate
    x_interpolate = get_neighbors(distance=distance_x, nearest_coord=nearest_x, spatial_res=spatial_res,
                                  spatial_sigfigs=spatial_sigfigs)  # figure out if nearest x coordinate is +/- of real x coordinate
    y_interpolate = get_neighbors(distance=distance_y, nearest_coord=nearest_y, spatial_res=spatial_res,
                                  spatial_sigfigs=spatial_sigfigs)  # figure out if nearest y coordinate is +/- of real y coordinate
    z_interpolate = get_neighbors(distance=distance_z, nearest_coord=nearest_z, spatial_res=spatial_res,
                                  spatial_sigfigs=spatial_sigfigs)  # figure out if nearest z coordinate is +/- of real z coordinate
    if x_interpolate < nearest_x:
        x_plus_coord = nearest_x
        x_minus_coord = x_interpolate
    elif x_interpolate > nearest_x:
        x_plus_coord = x_interpolate
        x_minus_coord = nearest_x
    elif x_interpolate <= 0.0:
        x_plus_coord = nearest_x
        x_minus_coord = max_x
    elif x_interpolate >= max_x:
        x_plus_coord = round(0.0, spatial_sigfigs)
        x_minus_coord = nearest_x
    else:
        x_plus_coord = nearest_x
        x_minus_coord = nearest_x

    if y_interpolate < nearest_y:
        y_plus_coord = nearest_y
        y_minus_coord = y_interpolate
    elif y_interpolate > nearest_y:
        y_plus_coord = y_interpolate
        y_minus_coord = nearest_y
    elif y_interpolate <= 0.0:
        y_plus_coord = nearest_y
        y_minus_coord = max_y
    elif y_interpolate >= max_y:
        y_plus_coord = round(0.0, spatial_sigfigs)
        y_minus_coord = nearest_y
    else:
        y_plus_coord = nearest_y
        y_minus_coord = nearest_y

    if z_interpolate < nearest_z:
        z_plus_coord = nearest_z
        z_minus_coord = z_interpolate
    elif z_interpolate > nearest_z:
        z_plus_coord = z_interpolate
        z_minus_coord = nearest_z
    elif z_interpolate <= 0.0:
        z_plus_coord = nearest_z
        z_minus_coord = max_z
    elif z_interpolate >= max_z:
        z_plus_coord = round(0.0, spatial_sigfigs)
        z_minus_coord = nearest_z
    else:
        z_plus_coord = nearest_z
        z_minus_coord = nearest_z

    possible_verteces = [
        (x_minus_coord, y_plus_coord, z_plus_coord),
        (x_minus_coord, y_plus_coord, z_minus_coord),
        (x_minus_coord, y_minus_coord, z_plus_coord),
        (x_minus_coord, y_minus_coord, z_minus_coord),
        (x_plus_coord, y_plus_coord, z_plus_coord),
        (x_plus_coord, y_plus_coord, z_minus_coord),
        (x_plus_coord, y_minus_coord, z_plus_coord),
        (x_plus_coord, y_minus_coord, z_minus_coord),
    ]

    cell_verteces = sorted(set(possible_verteces))

    cell_indices = []
    for i in cell_verteces:
        cell_indices.append(predict_index(coord=i, max_x=max_x, max_y=max_y,
                                           max_z=max_z, spatial_res=spatial_res, verbose=verbose))

    return nearest_index, cell_indices, cell_verteces

def update_position(coord, velocity, delta_time, max_z):
    """
    Accepts 3D coord, 3D velcoity vectors, and the timestep.  Returns new coordinate position.
    :param coord:
    :param velocity:
    :param delta_time:
    :return:
    """
    x, y, z = coord[0], coord[1], coord[2]
    x_vel, y_vel, z_vel = velocity[0], velocity[1], velocity[2]
    change_x, change_y, change_z = (x_vel * delta_time), (y_vel * delta_time), (z_vel * delta_time)
    new_x, new_y, new_z = (x + change_x), (y + change_y), (z + change_z)
    if new_z < max_z:
        new_coord = (new_x, new_y, new_z)
    else:
        new_coord = (new_x, new_y, max_z)
    return new_coord


def object_actions(mesh_df, objects_df, spatial_res, spatial_sigfigs, evolution_time, initial_time, max_x, max_y, max_z,
                   x_plus, x_minus, y_plus, y_minus, z_plus, z_minus, delta_time, matrix_densities, matrix_viscosities,
                   lower_model, upper_model, verbose=True):
    
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
    for object_index, object_object in enumerate(object_objects):
        object_id = object_object_ids[object_index]  # get the current object id
        coord = object_coords[object_index]  # get the current object coordinate
        #  interpolate to find the cell in which the object exists
        #  calculate the object's settling velocity
        if evolution_time == initial_time:
            initial_cell = interpolate_cell(coord=coord, spatial_res=spatial_res,
                                                     spatial_sigfigs=spatial_sigfigs, max_x=max_x,
                                                     max_y=max_y, max_z=max_z, x_plus=x_plus, x_minus=x_minus,
                                                     y_plus=y_plus, y_minus=y_minus, z_plus=z_plus, z_minus=z_minus,
                                                     verbose=verbose)
            nearest_indices[object_index] = initial_cell[0]
            cell_indices[object_index] = initial_cell[1]
            cell_vertices[object_index] = initial_cell[2]
        velocity = settling_modes.stokes_terminal(density=object_densities[object_index],
                                                  density_matrix=matrix_densities[nearest_indices[object_index]],
                                                  drag_coeff=2, radius=object_radii[object_index],
                                                  viscosity_matrix=matrix_viscosities[nearest_indices[object_index]],
                                                  current_coord=coord, lower_model=lower_model, upper_model=upper_model)
        #  get the object's new coordinates based on the object's velocity
        updated_coords = update_position(coord=coord, velocity=velocity, delta_time=delta_time,
                                                  max_z=max_z)
        cell = interpolate_cell(coord=updated_coords, spatial_res=spatial_res,
                                         spatial_sigfigs=spatial_sigfigs, max_x=max_x,
                                         max_y=max_y, max_z=max_z, x_plus=x_plus, x_minus=x_minus,
                                         y_plus=y_plus, y_minus=y_minus, z_plus=z_plus, z_minus=z_minus,
                                         verbose=verbose)
        console.event("{} ({}) will travel from {} to {} (velocity: {})".format(
            object_object, object_id, coord, updated_coords, velocity), verbose=verbose)
        #  update the dataframes with the new data
        object_velocities[object_index] = velocity
        object_coords[object_index] = updated_coords
        nearest_indices[object_index] = cell[0]
        cell_indices[object_index] = cell[1]
        cell_vertices[object_index] = cell[2]
    objects_df['coords'] = object_coords
    objects_df['velocity'] = object_velocities
    objects_df['cell_vertices'] = cell_vertices
    objects_df['cell_indices'] = cell_indices
    objects_df['nearest_index'] = nearest_indices
    return object_coords, nearest_indices, cell_indices
