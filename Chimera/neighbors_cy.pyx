

def get_neighbor_index(coord, array):
    index = array.index(coord)
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
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array)

        elif x_coord == min_x:
            x_plus = round(x_coord + spatial_res, spatial_sigfigs)
            x_minus = round(float(max_x), spatial_sigfigs)
            x_plus_coord = (x_plus, y_coord, z_coord)
            x_minus_coord = (x_minus, y_coord, z_coord)
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array)

        elif x_coord == max_x:
            x_plus = round(float(min_x), spatial_sigfigs)
            x_minus = round(x_coord - spatial_res, spatial_sigfigs)
            x_plus_coord = (x_plus, y_coord, z_coord)
            x_minus_coord = (x_minus, y_coord, z_coord)
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array)

        if min_y < y_coord < max_y:
            y_plus = round(y_coord + spatial_res, spatial_sigfigs)
            y_minus = round(y_coord - spatial_res, spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus, z_coord)
            y_minus_coord = (x_coord, y_minus, z_coord)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array)

        elif y_coord == min_y:
            y_plus = round(y_coord + spatial_res, spatial_sigfigs)
            y_minus = round(float(max_y), spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus, z_coord)
            y_minus_coord = (x_coord, y_minus, z_coord)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array)

        elif y_coord == max_y:
            y_plus = round(float(min_y), spatial_sigfigs)
            y_minus = round(y_coord - spatial_res, spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus, z_coord)
            y_minus_coord = (x_coord, y_minus, z_coord)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array)

        if min_z < z_coord < max_z:
            z_plus = round(z_coord + spatial_res, spatial_sigfigs)
            z_minus = round(z_coord - spatial_res, spatial_sigfigs)
            z_plus_coord = (x_coord, y_coord, z_plus)
            z_minus_coord = (x_coord, y_coord, z_minus)
            z_plus_index = get_neighbor_index(coord=z_plus_coord, array=array)
            z_minus_index = get_neighbor_index(coord=z_minus_coord, array=array)

        elif z_coord == min_z:
            z_plus = round(z_coord + spatial_res, spatial_sigfigs)
            z_minus = round(float(max_z), spatial_sigfigs)
            z_plus_coord = (x_coord, y_coord, z_plus)
            z_minus_coord = (x_coord, y_coord, z_minus)
            z_plus_index = get_neighbor_index(coord=z_plus_coord, array=array)
            z_minus_index = get_neighbor_index(coord=z_minus_coord, array=array)

        elif z_coord == max_z:
            z_plus = round(float(min_z), spatial_sigfigs)
            z_minus = round(z_coord - spatial_res, spatial_sigfigs)
            z_plus_coord = (x_coord, y_coord, z_plus)
            z_minus_coord = (x_coord, y_coord, z_minus)
            z_plus_index = get_neighbor_index(coord=z_plus_coord, array=array)
            z_minus_index = get_neighbor_index(coord=z_minus_coord, array=array)

        return x_plus_index, x_minus_index, y_plus_index, y_minus_index, z_plus_index, z_minus_index
    else:
        if min_x < x_coord < max_x:
            x_plus = round(x_coord + spatial_res, spatial_sigfigs)
            x_minus = round(x_coord - spatial_res, spatial_sigfigs)
            x_plus_coord = (x_plus, y_coord)
            x_minus_coord = (x_minus, y_coord)
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array)

        elif x_coord == min_x:
            x_plus = round(x_coord + spatial_res, spatial_sigfigs)
            x_minus = round(float(max_x), spatial_sigfigs)
            x_plus_coord = (x_plus, y_coord)
            x_minus_coord = (x_minus, y_coord)
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array)

        elif x_coord == max_x:
            x_plus = round(float(min_x), spatial_sigfigs)
            x_minus = round(x_coord - spatial_res, spatial_sigfigs)
            x_plus_coord = (x_plus, y_coord)
            x_minus_coord = (x_minus, y_coord)
            x_plus_index = get_neighbor_index(coord=x_plus_coord, array=array)
            x_minus_index = get_neighbor_index(coord=x_minus_coord, array=array)

        if min_y < y_coord < max_y:
            y_plus = round(y_coord + spatial_res, spatial_sigfigs)
            y_minus = round(y_coord - spatial_res, spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus)
            y_minus_coord = (x_coord, y_minus)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array)

        elif y_coord == min_y:
            y_plus = round(y_coord + spatial_res, spatial_sigfigs)
            y_minus = round(float(max_y), spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus)
            y_minus_coord = (x_coord, y_minus)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array)

        elif y_coord == max_y:
            y_plus = round(float(min_y), spatial_sigfigs)
            y_minus = round(y_coord - spatial_res, spatial_sigfigs)
            y_plus_coord = (x_coord, y_plus)
            y_minus_coord = (x_coord, y_minus)
            y_plus_index = get_neighbor_index(coord=y_plus_coord, array=array)
            y_minus_index = get_neighbor_index(coord=y_minus_coord, array=array)

        return x_plus_index, x_minus_index, y_plus_index, y_minus_index



