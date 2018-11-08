import os
import pandas as pd
import matplotlib as mpl; mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from ast import literal_eval


def interpolateCell(coord, spatial_sigfigs, spatial_res, max_x, max_y, max_z):
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
    nearest_x, nearest_y, nearest_z = arbitrary_round(
        coord=coord,
        spatial_res=spatial_res,
        spatial_sigfigs=spatial_sigfigs
    )  # find the nearest defined coordinates in the mesh
    nearest_coord = (nearest_x, nearest_y, nearest_z)  # the nearest coordinate point (tuple)
    distance_x, distance_y, distance_z = (x - nearest_x), (y - nearest_y), (z - nearest_z)  # find the distance of the real coordinate to the nearest mesh coordinate
    x_interpolate = get_neighbors(
        distance=distance_x,
        nearest_coord=nearest_x,
        spatial_res=spatial_res,
        spatial_sigfigs=spatial_sigfigs
    )  # figure out if nearest x coordinate is +/- of real x coordinate
    y_interpolate = get_neighbors(
        distance=distance_y,
        nearest_coord=nearest_y,
        spatial_res=spatial_res,
        spatial_sigfigs=spatial_sigfigs
    )  # figure out if nearest y coordinate is +/- of real y coordinate
    z_interpolate = get_neighbors(
        distance=distance_z,
        nearest_coord=nearest_z,
        spatial_res=spatial_res,
        spatial_sigfigs=spatial_sigfigs
    )  # figure out if nearest z coordinate is +/- of real z coordinate
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

    return cell_verteces


def parseMesh(const_xval, const_yval, min_zval, max_zval, spatial_res, spatial_sigfigs, max_x, max_y, max_z):
    df = pd.read_csv("mesh_.csv")
    for row in df.index:
        coords_tuple = literal_eval(df['coords'][row])
        if min_zval <= coords_tuple[2] <= max_zval:
            cell_vertices = interpolateCell(coord=coords_tuple, spatial_sigfigs=spatial_sigfigs, spatial_res=spatial_res,
                            max_x=max_x, max_y=max_y, max_z=max_z)





if __name__ == "__main__":
    const_xval = 1.2
    const_yval = 2.3
    min_zval = 1.0
    max_zval = 10.0
    spatial_res = 0.5
    spatial_sigfigs = 2
    max_x = 5.0
    max_y = 5.0
    max_z = 15.0
    parseMesh(const_xval=const_xval, const_yval=const_yval, min_zval=min_zval, max_zval=max_zval,
              spatial_res=spatial_res, spatial_sigfigs=spatial_sigfigs, max_x=max_x, max_y=max_y, max_z=max_z)



    # print("Please enter the x-value to hold constant.")
    # const_xval = input(">>> ")
    # print("Please enter the y-value to hold constant.")
    # const_yval = input(">>> ")
    # print("Please enter the minimum z-value.")
    # min_zval = input(">>> ")
    # print("Please enter the maximum z-value.")
    # max_zval = input(">>> ")
    # print("Please enter the spatial resolution of the model.")
    # spatial_res = input(">>> ")
    #
    # print("Would you like to strip 1 file (one) or a series of files (series)?")
    # strip_option = input(">>> ").lower()

    # if strip_option == "one":
        # parseMesh(const_xval=const_xval, const_yval=const_yval, min_zval=min_zval,
        #           max_zval=max_zval, spatial_res=spatial_res)