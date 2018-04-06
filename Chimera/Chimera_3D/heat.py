import multiprocessing as mp
from math import pi
from . import backends, dynamics

def conduction(coords, len_coords, x_plus_indices, x_minus_indices, y_plus_indices, y_minus_indices, z_plus_indices,
               z_minus_indices, temperatures, object_ids, spatial_res, conductivities, delta_time, mesh_indices,
               num_workers, multiprocess=False):

    if multiprocess is False:

        update_temps = [0 for _ in range(0, len_coords)]  # all of the updated temperatures due to conduction
        dT_dt_list = [0 for _ in range(0, len_coords)]

        for index, coord in enumerate(coords):
            temp_point = temperatures[index]  # temperature of the specified coordinate position
            if 'C' not in object_ids[index]:  # makes sure that point z is not a boundary layer
                x_plus_index = x_plus_indices[index]  # index of the x+ coordinate position
                x_minus_index = x_minus_indices[index]  # index of the x- coordinate position
                y_plus_index = y_plus_indices[index]  # index of the y+ coordinate position
                y_minus_index = y_minus_indices[index]  # index of the y- coordinate position
                z_plus_index = z_plus_indices[index]  # index of the z+ coordinate position
                z_minus_index = z_minus_indices[index]  # index of the z- coordinate position
                temp_x_plus = temperatures[x_plus_index]  # temperature of the x+ coordinate position
                temp_x_minus = temperatures[x_minus_index]  # temperature of the x- coordinate position
                temp_y_plus = temperatures[y_plus_index]  # temperature of the y+ coordinate position
                temp_y_minus = temperatures[y_minus_index]  # temperature of the y- coordinate position
                temp_z_plus = temperatures[z_plus_index]  # temperature of the z+ coordinate position
                temp_z_minus = temperatures[z_minus_index]  # temperature of the z- coordinate position
                k = conductivities[index]  # conductivities of the material at position of coordinate z

                # central difference laplacian is ((f_(x-1) - (2*f_(x)) + f_(x+1)) / (delta_x)^2)
                # central difference laplacian is as follows for each vector component
                x_temp_laplacian = ((temp_x_plus - (2 * temp_point) + temp_x_minus) / ((spatial_res) ** 2))
                y_temp_laplacian = ((temp_y_plus - (2 * temp_point) + temp_y_minus) / ((spatial_res) ** 2))
                z_temp_laplacian = ((temp_z_plus - (2 * temp_point) + temp_z_minus) / ((spatial_res) ** 2))
                temp_laplacian = x_temp_laplacian + y_temp_laplacian + z_temp_laplacian

                # change in temperature with respect to time, dT/dt = -k * laplacian(T)
                dT_dt = (k) * temp_laplacian  # the central finite difference heat equation
                dT_dt_list[index] = dT_dt

                dT = dT_dt * delta_time  # the change in temperature with respect to the finite timestep
                new_T = temp_point + dT  # adds dT to the original temperature
                update_temps[index] = new_T  # adds the new temperature to the updated temperature list
            else:  # if it is a boundary layer, it is a fixed temperature
                update_temps[index] = temp_point
                dT_dt_list[index] = 0.0

        return update_temps, dT_dt_list

    else:
        update_temps, dT_dt_list = multiprocess_conduction_manager(coords=coords, len_coords=len_coords,
                                                                   x_plus_indices=x_plus_indices,
                                                                   x_minus_indices=x_minus_indices,
                                                                   y_plus_indices=y_plus_indices,
                                                                   y_minus_indices=y_minus_indices,
                                                                   z_plus_indices=z_plus_indices,
                                                                   z_minus_indices=z_minus_indices,
                                                                   temperatures=temperatures,
                                                                   conductivities=conductivities,
                                                                   spatial_res=spatial_res, delta_time=delta_time,
                                                                   object_ids=object_ids, mesh_indices=mesh_indices,
                                                                   num_workers=num_workers,
                                                                   multiprocess=multiprocess)
        return update_temps, dT_dt_list

def multiprocess_conduction(coords, len_coords, x_plus_indices, x_minus_indices, y_plus_indices, y_minus_indices, z_plus_indices,
               z_minus_indices, temperatures, object_ids, spatial_res, conductivities, delta_time, mesh_indices,
               num_workers, update_temps, dT_dt_list, multiprocess=True):

    for index, coord in enumerate(coords):
        temp_point = temperatures[index]  # temperature of the specified coordinate position
        if 'C' not in object_ids[index]:  # makes sure that point z is not a boundary layer
            x_plus_index = x_plus_indices[index]  # index of the x+ coordinate position
            x_minus_index = x_minus_indices[index]  # index of the x- coordinate position
            y_plus_index = y_plus_indices[index]  # index of the y+ coordinate position
            y_minus_index = y_minus_indices[index]  # index of the y- coordinate position
            z_plus_index = z_plus_indices[index]  # index of the z+ coordinate position
            z_minus_index = z_minus_indices[index]  # index of the z- coordinate position
            temp_x_plus = temperatures[x_plus_index]  # temperature of the x+ coordinate position
            temp_x_minus = temperatures[x_minus_index]  # temperature of the x- coordinate position
            temp_y_plus = temperatures[y_plus_index]  # temperature of the y+ coordinate position
            temp_y_minus = temperatures[y_minus_index]  # temperature of the y- coordinate position
            temp_z_plus = temperatures[z_plus_index]  # temperature of the z+ coordinate position
            temp_z_minus = temperatures[z_minus_index]  # temperature of the z- coordinate position
            k = conductivities[index]  # conductivities of the material at position of coordinate z

            # central difference laplacian is ((f_(x-1) - (2*f_(x)) + f_(x+1)) / (delta_x)^2)
            # central difference laplacian is as follows for each vector component
            x_temp_laplacian = ((temp_x_plus - (2 * temp_point) + temp_x_minus) / ((spatial_res) ** 2))
            y_temp_laplacian = ((temp_y_plus - (2 * temp_point) + temp_y_minus) / ((spatial_res) ** 2))
            z_temp_laplacian = ((temp_z_plus - (2 * temp_point) + temp_z_minus) / ((spatial_res) ** 2))
            temp_laplacian = x_temp_laplacian + y_temp_laplacian + z_temp_laplacian

            # change in temperature with respect to time, dT/dt = -k * laplacian(T)
            dT_dt = (k) * temp_laplacian  # the central finite difference heat equation
            dT_dt_list[index] = dT_dt

            dT = dT_dt * delta_time  # the change in temperature with respect to the finite timestep
            new_T = temp_point + dT  # adds dT to the original temperature
            update_temps[index] = new_T  # adds the new temperature to the updated temperature list
        else:  # if it is a boundary layer, it is a fixed temperature
            update_temps[index] = temp_point
            dT_dt_list[index] = 0.0

    return update_temps, dT_dt_list


def multiprocess_conduction_manager(coords, len_coords, x_plus_indices, x_minus_indices, y_plus_indices, y_minus_indices, z_plus_indices,
               z_minus_indices, temperatures, object_ids, spatial_res, conductivities, delta_time, mesh_indices,
               num_workers, multiprocess=True):


    update_temps = mp.Array('d', len_coords)  # all of the updated temperatures due to conduction
    dT_dt_list = mp.Array('d', len_coords)

    jobs = []

    index_chunks = backends.chunk_array(array=mesh_indices, num_chunks=num_workers, len_array=len_coords)
    for process_id, process in enumerate(index_chunks):
        start_index, stop_index = process[0], process[-1]
        chunk_coords = coords[start_index:stop_index]
        chunk_x_plus_indices = x_plus_indices[start_index:stop_index]
        chunk_x_minus_indices = x_minus_indices[start_index:stop_index]
        chunk_y_plus_indices = y_plus_indices[start_index:stop_index]
        chunk_y_minus_indices = y_minus_indices[start_index:stop_index]
        chunk_z_plus_indices = z_plus_indices[start_index:stop_index]
        chunk_z_minus_indices = z_minus_indices[start_index:stop_index]
        chunk_object_ids = object_ids[start_index:stop_index]
        p = mp.Process(name=str(process_id), target=multiprocess_conduction, args=(chunk_coords, len_coords,
                                                                                   chunk_x_plus_indices,
                                                                                   chunk_x_minus_indices,
                                                                                   chunk_y_plus_indices,
                                                                                   chunk_y_minus_indices,
                                                                                   chunk_z_plus_indices,
                                                                                   chunk_z_minus_indices,
                                                                                   temperatures,
                                                                                   chunk_object_ids, spatial_res,
                                                                                   conductivities,
                                                                                   delta_time, mesh_indices,
                                                                                   num_workers, update_temps,
                                                                                   dT_dt_list, multiprocess))
        jobs.append(p)
        # p.start()
    for job in jobs:
        job.start()
        job.join()

    return update_temps, dT_dt_list

def object_conduction(object_temperatures, object_index, object_k, mesh_temperatures, nearest_index, farthest_index,
                      spatial_res, delta_time, distances, directional_vertices, vertex_distances):

    nearest_temp = mesh_temperatures[nearest_index]
    farthest_temp = mesh_temperatures[farthest_index]
    print("\n", nearest_temp, object_temperatures[object_index])
    t_laplacian = (((nearest_temp) - (2 * object_temperatures[object_index]) + farthest_temp) / ((spatial_res) ** 2))
    dT_dt = object_k * t_laplacian
    dT = dT_dt * delta_time
    if dT > 0:
        mesh_temperatures[nearest_index] -= dT
        object_temperatures[object_index] += dT
        print("\nMESH WILL LOSE {}, OBJECT WILL GAIN {}\nMESH T: {}, OBJECT T: {}".format(-dT, dT, mesh_temperatures[nearest_index], mesh_temperatures[object_index]))
    elif dT < 0:
        mesh_temperatures[nearest_index] -= dT
        object_temperatures[object_index] += dT
        print("\nMESH WILL GAIN {}, OBJECT WILL LOSE {}\nMESH T: {}, OBJECT T: {}".format(-dT, dT, mesh_temperatures[nearest_index], mesh_temperatures[object_index]))

    else:
        pass
    return dT, dT_dt


def viscous_dissipation(drag_coeff, matrix_density, object_radius, object_density, object_velocity, delta_time, cp,
                        distance_travelled):
    """
    The heat in Kelvin produced by viscous dissipation.
    :param drag_coeff:
    :param matrix_density:
    :param object_radius:
    :param object_velocity:
    :param object_mass:
    :param delta_time:
    :param cp:
    :return:
    """
    object_volume = (4/3) * pi * (object_radius**3)
    object_mass = object_density * object_volume
    grav_accel = 9.81
    # the drag force on the object, which will be converted to heat via viscous dissipation
    fd = dynamics.drag_force(drag_coeff=drag_coeff, matrix_density=matrix_density, object_radius=object_radius,
                             distance_travelled=distance_travelled, delta_time=delta_time)
    fb = dynamics.buoyant_force(matrix_density=matrix_density, grav_accel=grav_accel, object_radius=object_radius)
    fg = dynamics.grav_force(object_radius=object_radius, object_density=object_density, grav_accel=grav_accel)
    # for simplicity, assume the drag force is conservative within the time interval
    # use this conservative drag force to find the work done by it
    w = dynamics.work_conservative(delta_time=delta_time, force=fd, velocity=object_velocity, distance_travelled=distance_travelled)
    # convert the work to degrees kelvin
    k = w / (cp * object_mass)
    return k
