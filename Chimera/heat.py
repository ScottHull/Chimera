from . import backends

def conduction_linear(z, z_plus_indeces, z_minus_indeces, temperatures, object_ids, spatial_res,
                      conductivity, delta_time):
    """
    Heat transfer through conduction alone.
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
    dT_dt_list = []
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
            dT_dt_list.append(dT_dt)

            dT = dT_dt * delta_time  # the change in temperature with respect to the finite timestep
            new_T = temp_z + dT  # adds dT to the original temperature
            update_temps.append(new_T)  # adds the new temperature to the updated temperature list
        else:
            update_temps.append(temp_z)
            dT_dt_list.append(0.0)

    return update_temps, dT_dt_list

def nusselt(layers, mesh, spatial_res):
    """
    The Nusselt Number provides the ratio of convective to conductive heat transfer.
    Nu = (convective flux / conductive flux) = (h/(k/L)) = ((h*L)/k)
    Where:
    h = convective heat transfer coefficient of the flow
    L = characteristic length
    k = thermal conductivity of the fluid
    :param layers:
    :return:
    """
    objects = layers['object'].tolist()
    z_min_coords = layers['min_z'].tolist()
    z_max_coords = layers['max_z'].tolist()
    nusselt_nos = layers['nusselt'].tolist()
    temperatures = mesh['temperature'].tolist()
    conductivities = mesh['conductivity'].tolist()
    dT_dts = mesh['dT_dt'].tolist()

    for index, object in enumerate(objects):  # loop through layers to determine
        z_min_coord = z_min_coords[index]
        z_min_index = backends.predict_linear_index(z=z_min_coord, spatial_res=spatial_res)
        z_max_coord = z_max_coords[index]
        z_max_index = backends.predict_linear_index(z=z_max_coord, spatial_res=spatial_res)
        z_min_temp = temperatures[z_min_index]
        z_max_temp = temperatures[z_max_index]
        if z_max_temp - z_min_temp != 0.0:
            z_max_dT_dt = dT_dts[z_max_index]
            z_min_dT_dt = dT_dts[z_min_index]
            # z_max_conductivity = conductivities[z_max_index]
            z_min_conductivity = conductivities[z_min_index]
            # z_max_h = backends.convect_heat_transfer_coeff(heat_flux=z_max_dT_dt,
            #                                                temp_bottom=z_min_temp, temp_top=z_max_temp)
            z_min_h = backends.convect_heat_transfer_coeff(heat_flux=z_min_dT_dt,
                                                           temp_bottom=z_max_temp, temp_top=z_min_temp)
            characteristic_length = z_max_coord - z_min_coord
            nusselt = (z_min_h * characteristic_length) / z_min_conductivity
            nusselt_nos[index] = nusselt
        else:
            nusselt_nos[index] = 0.0

    layers['nusselt'] = nusselt_nos

def rayleigh(density, gravity, thermal_expansivity, thermal_conductivity, temp_upper_layer, temp_lower_layer,
             length, viscosity, material_specific_heat):
    # Ra = (rho * g * alpha * delta_T * d^3) / (mu * k)
    delta_T = temp_lower_layer - temp_upper_layer
    diffusivity = backends.calculate_thermal_diffusivity(thermal_conductivity=thermal_conductivity, density=density,
                                                         specific_heat_capacity=material_specific_heat)
    rayleigh = (density * gravity * thermal_expansivity * delta_T * (length**3)) / (viscosity * diffusivity)
