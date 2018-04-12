from Chimera.Chimera_3D import console, dynamics, heat
import numpy as np
from math import sqrt

def set_columns(mesh_df, object_df, coords_range):

    mesh_columns = ["object", "object_id", "xplus_index", "xminus_index", "yplus_index", "yminus_index",
                    "zplus_index", "zminus_index", "temperature", "dT_dt", "conductivity", "viscosity", "density",
                    "diffusivity"]
    object_columns = ["object", "object_id", "coords", "temperature", "radius", "conductivity", "density", "velocity",
                      "cell_indices", "cell_vertices", "nearest_index", "drag_coeff", "cp", "diffusivity"]
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

    if object_type == 'matrix':
        object_id = random_gen(object_identifier='B', id_val=id_val)
        # while object_id in object_ids:
        #     object_id = random_gen(object_identifier='B')
        return object_id
    elif object_type == 'object':
        object_id = random_gen(object_identifier='A', id_val=id_val)
        return object_id
    elif object_type == 'boundary':
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
        return None

def override_timestep(timestep, conductivities, spatial_res, spatial_sigfigs, diffusivities, verbose):
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
        min_diffusivity = min(diffusivities[diffusivities != 0.0])
        delta_time = (spatial_res**2) / (2 * (min_diffusivity))
        delta_time = round(delta_time, spatial_sigfigs)
        return delta_time
    else:
        min_diffusivity = min(diffusivities[diffusivities != 0.0])
        delta_time_test = (spatial_res ** 2) / (2 * (min_diffusivity))
        delta_time_test = round(delta_time_test, spatial_sigfigs)
        delta_time = timestep
        if delta_time > delta_time_test:
            console.error("Warning!  Your selected timestep is greater than the threshold!", verbose=verbose)
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

    def determine_directional_indices(x, y, z, vertex, directional_dict, vertex_index):
        vertex_x, vertex_y, vertex_z = vertex[0], vertex[1], vertex[2]
        if vertex_x > x:
            directional_dict['x+'].append(vertex_index)
        elif vertex_x < x:
            directional_dict['x-'].append(vertex_index)
        else:
            directional_dict['x+'].append(vertex_index)
            directional_dict['x-'].append(vertex_index)
        if vertex_y > y:
            directional_dict['y+'].append(vertex_index)
        elif vertex_y < y:
            directional_dict['y-'].append(vertex_index)
        else:
            directional_dict['y+'].append(vertex_index)
            directional_dict['y-'].append(vertex_index)
        if vertex_z > z:
            directional_dict['z+'].append(vertex_index)
        elif vertex_z < z:
            directional_dict['z-'].append(vertex_index)
        else:
            directional_dict['z+'].append(vertex_index)
            directional_dict['z-'].append(vertex_index)
        return directional_dict

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
        x_minus_coord = round(nearest_x - spatial_res, spatial_sigfigs)

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
        y_minus_coord = round(nearest_y - spatial_res, spatial_sigfigs)

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
        z_minus_coord = round(nearest_z - spatial_res, spatial_sigfigs)

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

    total_distance = 0
    max_distance = 0
    farthest_coord = nearest_coord
    directional_vertices = {
        'x+': [],
        'x-': [],
        'y+': [],
        'y-': [],
        'z+': [],
        'z-': [],
    }

    vertex_distances = {
    }

    cell_indices = []
    for i in cell_verteces:
        vertex_index = predict_index(coord=i, max_x=max_x, max_y=max_y,
                      max_z=max_z, spatial_res=spatial_res, verbose=verbose)
        cell_indices.append(vertex_index)
        # distance = sqrt(((nearest_coord[0] - i[0])**2) + ((nearest_coord[1] - i[1])**2) + ((nearest_coord[2] - i[2])**2))
        distance = sqrt(((x - i[0])**2) + ((y - i[1])**2) + ((z - i[2])**2))
        vertex_distances.update({vertex_index: distance})
        total_distance += distance
        if distance > max_distance:
            max_distance = distance
            farthest_coord = i
        determine_directional_indices(x=x, y=y, z=z, directional_dict=directional_vertices, vertex=i, vertex_index=vertex_index)
    farthest_index = predict_index(coord=farthest_coord, max_x=max_x, max_y=max_y, max_z=max_z, spatial_res=spatial_res)

    return nearest_index, cell_indices, cell_verteces, farthest_index, directional_vertices, \
                                                            vertex_distances, total_distance

def update_position(coord, velocity, delta_time, max_z, spatial_res):
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
        distance_travelled = ((new_x - x), (new_y - y), (new_z - z))
    else:
        new_coord = (new_x, new_y, max_z)
        distance_travelled = (0, 0, 0)
    return new_coord, distance_travelled

def chunk_array(array, num_chunks, len_array):
    """
    Chunk a given list into a given number of smaller lists.
    :param array:
    :param num_chunks:
    :return:
    """
    # return [array[i:i + num_chunks] for i in range(0, len(array), num_chunks)]
    if len_array > num_chunks:
        if len_array % num_chunks == 0:
            chunks = list(zip(*[iter(array)] * (int(round(len_array / num_chunks)))))
            return chunks
        else:
            odd_num = (array[-1],)
            chunks = list(zip(*[iter(array)] * (int(round(len_array / num_chunks)))))
            chunks[-1] = chunks[-1] + odd_num
            return chunks
    else:
        return array