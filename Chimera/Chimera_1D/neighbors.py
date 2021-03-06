from . import backends


def get_neighbor_index(coord, array, max_x, max_y, max_z, verbose, spatial_res):
    # index = array.index(coord)
    # return index
    index = backends.predict_index(coord=coord, max_x=max_x, max_y=max_y, max_z=max_z, spatial_res=spatial_res,
                                   verbose=verbose)
    return index

def get_neighbors(verbose, coords, array, spatial_res, spatial_sigfigs, max_x, max_y,
                  max_z=False, min_x=0.0, min_y=0.0, min_z=0.0):

    if max_z is not False:
        x_coord, y_coord, z_coord = coords[0], coords[1], coords[2]
    else:
        x_coord, y_coord = coords[0], coords[1]
    if max_z is not False:
        if min_x < x_coord < max_x:
            x_plus = round(x_coord + spatial_res, spatial_sigfigs)
            x_minus = round(x_coord - spatial_res, spatial_sigfigs)
            x_plus_coord = (x_plus, y_coord, z_coord)
            x_minus_coord = (x_minus, y_coord, z_coord)
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        elif x_coord == min_x:
            x_plus = round(x_coord + spatial_res, spatial_sigfigs)
            x_minus = round(float(max_x), spatial_sigfigs)
            x_plus_coord = (x_plus, y_coord, z_coord)
            x_minus_coord = (x_minus, y_coord, z_coord)
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        elif x_coord == max_x:
            x_plus = round(float(min_x), spatial_sigfigs)
            x_minus = round(x_coord - spatial_res, spatial_sigfigs)
            x_plus_coord = (x_plus, y_coord, z_coord)
            x_minus_coord = (x_minus, y_coord, z_coord)
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        if min_y < y_coord < max_y:
            y_plus = round(y_coord + spatial_res, spatial_sigfigs)
            y_minus = round(y_coord - spatial_res, spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus, z_coord)
            y_minus_coord = (x_coord, y_minus, z_coord)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        elif y_coord == min_y:
            y_plus = round(y_coord + spatial_res, spatial_sigfigs)
            y_minus = round(float(max_y), spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus, z_coord)
            y_minus_coord = (x_coord, y_minus, z_coord)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        elif y_coord == max_y:
            y_plus = round(float(min_y), spatial_sigfigs)
            y_minus = round(y_coord - spatial_res, spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus, z_coord)
            y_minus_coord = (x_coord, y_minus, z_coord)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        if min_z < z_coord < max_z:
            z_plus = round(z_coord + spatial_res, spatial_sigfigs)
            z_minus = round(z_coord - spatial_res, spatial_sigfigs)
            z_plus_coord = (x_coord, y_coord, z_plus)
            z_minus_coord = (x_coord, y_coord, z_minus)
            z_plus_index = get_neighbor_index(coord=z_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            z_minus_index = get_neighbor_index(coord=z_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        elif z_coord == min_z:
            z_plus = round(z_coord + spatial_res, spatial_sigfigs)
            z_minus = round(float(max_z), spatial_sigfigs)
            z_plus_coord = (x_coord, y_coord, z_plus)
            z_minus_coord = (x_coord, y_coord, z_minus)
            z_plus_index = get_neighbor_index(coord=z_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            z_minus_index = get_neighbor_index(coord=z_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        elif z_coord == max_z:
            z_plus = round(float(min_z), spatial_sigfigs)
            z_minus = round(z_coord - spatial_res, spatial_sigfigs)
            z_plus_coord = (x_coord, y_coord, z_plus)
            z_minus_coord = (x_coord, y_coord, z_minus)
            z_plus_index = get_neighbor_index(coord=z_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            z_minus_index = get_neighbor_index(coord=z_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        return x_plus_index, x_minus_index, y_plus_index, y_minus_index, z_plus_index, z_minus_index
    else:
        if min_x < x_coord < max_x:
            x_plus = round(x_coord + spatial_res, spatial_sigfigs)
            x_minus = round(x_coord - spatial_res, spatial_sigfigs)
            x_plus_coord = (x_plus, y_coord)
            x_minus_coord = (x_minus, y_coord)
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        elif x_coord == min_x:
            x_plus = round(x_coord + spatial_res, spatial_sigfigs)
            x_minus = round(float(max_x), spatial_sigfigs)
            x_plus_coord = (x_plus, y_coord)
            x_minus_coord = (x_minus, y_coord)
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        elif x_coord == max_x:
            x_plus = round(float(min_x), spatial_sigfigs)
            x_minus = round(x_coord - spatial_res, spatial_sigfigs)
            x_plus_coord = (x_plus, y_coord)
            x_minus_coord = (x_minus, y_coord)
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        if min_y < y_coord < max_y:
            y_plus = round(y_coord + spatial_res, spatial_sigfigs)
            y_minus = round(y_coord - spatial_res, spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus)
            y_minus_coord = (x_coord, y_minus)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        elif y_coord == min_y:
            y_plus = round(y_coord + spatial_res, spatial_sigfigs)
            y_minus = round(float(max_y), spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus)
            y_minus_coord = (x_coord, y_minus)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        elif y_coord == max_y:
            y_plus = round(float(min_y), spatial_sigfigs)
            y_minus = round(y_coord - spatial_res, spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus)
            y_minus_coord = (x_coord, y_minus)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array, max_x=max_x, max_y=max_y,
                                              max_z=max_z, spatial_res=spatial_res, verbose=verbose)

        return x_plus_index, x_minus_index, y_plus_index, y_minus_index


def get_linear_neighbor_index(coord, spatial_res, spatial_sigfigs):
    """
    Returns the index position of the specified coordinate in a linear model.
    :param coord:
    :param spatial_res:
    :param spatial_sigfigs:
    :return:
    """

    index = int(round((coord / spatial_res), spatial_sigfigs))
    return index


def get_linear_neighbors(verbose, coord, spatial_res, spatial_sigfigs,
                  max_z, min_z=0.0):
    """
    Fetch the two neighboring coordinate index positions of the specified coordinate.
    :param verbose:
    :param coord:
    :param spatial_res:
    :param spatial_sigfigs:
    :param max_z:
    :param min_z:
    :return:
    """

    if float(coord) == float(max_z):
        z_plus = round(min_z, spatial_sigfigs)
        z_minus = round(coord - spatial_res, spatial_sigfigs)
        z_plus_index = get_linear_neighbor_index(z_plus, spatial_res=spatial_res, spatial_sigfigs=spatial_sigfigs)
        z_minus_index = get_linear_neighbor_index(z_minus, spatial_res=spatial_res, spatial_sigfigs=spatial_sigfigs)
    elif float(coord) == float(min_z):
        z_plus = round(coord + spatial_res, spatial_sigfigs)
        z_minus = round(max_z, spatial_sigfigs)
        z_plus_index = get_linear_neighbor_index(z_plus, spatial_res=spatial_res, spatial_sigfigs=spatial_sigfigs)
        z_minus_index = get_linear_neighbor_index(z_minus, spatial_res=spatial_res, spatial_sigfigs=spatial_sigfigs)
    else:
        z_plus = round(coord + spatial_res, spatial_sigfigs)
        z_minus = round(coord - spatial_res, spatial_sigfigs)
        z_plus_index = get_linear_neighbor_index(z_plus, spatial_res=spatial_res, spatial_sigfigs=spatial_sigfigs)
        z_minus_index = get_linear_neighbor_index(z_minus, spatial_res=spatial_res, spatial_sigfigs=spatial_sigfigs)

    return z_plus_index, z_minus_index
