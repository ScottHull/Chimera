
def conduction(coords, x_plus_indeces, x_minus_indeces, y_plus_indeces, y_minus_indeces, z_plus_indeces,
               z_minus_indeces, temperatures, object_ids, spatial_res, conductivity, delta_time):

    update_temps = []  # all of the updated temperatures due to conduction
    dT_dt_list = []
    for index, coord in enumerate(coords):
        temp_point = temperatures[index]  # temperature of the specified coordinate position
        if 'C' not in object_ids[index]:  # makes sure that point z is not a boundary layer
            x_plus_index = x_plus_indeces[index]  # index of the x+ coordinate position
            x_minus_index = x_minus_indeces[index]  # index of the x- coordinate position
            y_plus_index = y_plus_indeces[index]  # index of the y+ coordinate position
            y_minus_index = y_minus_indeces[index]  # index of the y- coordinate position
            z_plus_index = z_plus_indeces[index]  # index of the z+ coordinate position
            z_minus_index = z_minus_indeces[index]  # index of the z- coordinate position
            temp_x_plus = temperatures[x_plus_index]  # temperature of the x+ coordinate position
            temp_x_minus = temperatures[x_minus_index]  # temperature of the x- coordinate position
            temp_y_plus = temperatures[y_plus_index]  # temperature of the y+ coordinate position
            temp_y_minus = temperatures[y_minus_index]  # temperature of the y- coordinate position
            temp_z_plus = temperatures[z_plus_index]  # temperature of the z+ coordinate position
            temp_z_minus = temperatures[z_minus_index]  # temperature of the z- coordinate position
            k = conductivity[index]  # conductivity of the material at position of coordinate z

            # central difference laplacian is ((f_(x-1) - (2*f_(x)) + f_(x+1)) / (delta_x)^2)
            # central difference laplacian is as follows for each vector component
            x_temp_laplacian = ((temp_x_plus - (2 * temp_point) + temp_x_minus) / ((spatial_res) ** 2))
            y_temp_laplacian = ((temp_y_plus - (2 * temp_point) + temp_y_minus) / ((spatial_res) ** 2))
            z_temp_laplacian = ((temp_z_plus - (2 * temp_point) + temp_z_minus) / ((spatial_res) ** 2))
            temp_laplacian = x_temp_laplacian + y_temp_laplacian + z_temp_laplacian

            # change in temperature with respect to time, dT/dt = -k * laplacian(T)
            dT_dt = (k) * temp_laplacian  # the central finite difference heat equation
            dT_dt_list.append(dT_dt)

            dT = dT_dt * delta_time  # the change in temperature with respect to the finite timestep
            new_T = temp_point + dT  # adds dT to the original temperature
            update_temps.append(new_T)  # adds the new temperature to the updated temperature list
        else:  # if it is a boundary layer, it is a fixed temperature
            update_temps.append(temp_point)
            dT_dt_list.append(0.0)

    return update_temps, dT_dt_list


def viscous_dissipation():
    pass

