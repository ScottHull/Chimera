from . import backends

def conduction_linear(z, z_plus_indeces, z_minus_indeces, temperatures, object_ids, spatial_res,
                      conductivity, delta_time):

    """
    Heat transfer through conduction.
    :param z:
    :param z_plus_indeces:
    :param z_minus_indeces:
    :param temperatures:
    :param spatial_res:
    :param conductivity:
    :param delta_time:
    :return:
    """

    update_temps = []  # all of the updated temperatures due to conduction
    for index, coord in enumerate(z):
        temp_z = temperatures[index]  # temperature of the z coordinate position
        if 'C' not in object_ids[index]:  # makes sure that point z is not a boundary layer
            z_plus_index = z_plus_indeces[index]  # index of the z+ coordinate position
            z_minus_index = z_minus_indeces[index]  # index of the z- coordinate position
            temp_z_plus = temperatures[z_plus_index]  # temperature of the z+ coordinate position
            temp_z_minus = temperatures[z_minus_index]  # temperature of the z- coordinate position
            k = conductivity[index]  # conductivity of the material at position of coordinate z

            # central difference laplacian is ((f_(x-1) - (2*f_(x)) + f_(x+1)) / (delta_x)^2)
            # central difference laplacian is as follows
            temp_laplacian = ((temp_z_plus - (2 * temp_z) + temp_z_minus) / ((spatial_res) ** 2))

            # change in temperature with respect to time, dT/dt = -k * laplacian(T)
            dT_dt = (k) * temp_laplacian  # the central finite difference heat equation

            dT = dT_dt * delta_time  # the change in temperature with respect to the finite timestep
            new_T = temp_z + dT  # adds dT to the original temperature
            update_temps.append(new_T)  # adds the new temperature to the updated temperature list
        else:
            update_temps.append(temp_z)

    return update_temps