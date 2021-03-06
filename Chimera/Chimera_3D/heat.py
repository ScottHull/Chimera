import multiprocessing as mp
from math import pi
import numpy as np
from . import backends, dynamics

def conduction(coords, len_coords, x_plus_indices, x_minus_indices, y_plus_indices, y_minus_indices, z_plus_indices,
               z_minus_indices, temperatures, object_ids, spatial_res, conductivities, delta_time, real_delta_time,
               mesh_indices, num_workers, multiprocess=False):

    if not multiprocess:

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
                dT_dt = k * temp_laplacian  # the central finite difference heat equation
                dT_dt_list[index] = dT_dt

                dT = dT_dt * (delta_time / real_delta_time)  # the change in temperature with respect to the finite normalized timestep
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

def object_conduction(object_temperatures, object_index, object_k, matrix_k, mesh_temperatures, nearest_index, farthest_index,
                      spatial_res, delta_time, directional_vertices, vertex_distances, total_distance, matrix_ids,
                      coords, matrix_diffusivities, spatial_sigfigs, verbose):

    box_spatial_res = (spatial_res / 2)
    box_delta_time = backends.override_timestep(spatial_res=box_spatial_res, conductivities=None,
                                                diffusivities=np.array([matrix_diffusivities[nearest_index]]),
                                                spatial_sigfigs=spatial_sigfigs, timestep=False, verbose=verbose)
    object_temp = object_temperatures[object_index]
    x_plus_temps = [mesh_temperatures[i] for i in directional_vertices["x+"]]
    x_minus_temps = [mesh_temperatures[i] for i in directional_vertices["x-"]]
    y_plus_temps = [mesh_temperatures[i] for i in directional_vertices["y+"]]
    y_minus_temps = [mesh_temperatures[i] for i in directional_vertices["y-"]]
    z_plus_temps = [mesh_temperatures[i] for i in directional_vertices["z+"]]
    z_minus_temps = [mesh_temperatures[i] for i in directional_vertices["z-"]]
    d2T_dx2 = ((sum(x_plus_temps) / len(x_plus_temps)) - (2 * object_temp) + (sum(x_minus_temps) / len(x_minus_temps))) / (
                (box_spatial_res ** 2))
    d2T_dy2 = ((sum(y_plus_temps) / len(y_plus_temps)) - (2 * object_temp) + (sum(y_minus_temps) / len(y_minus_temps))) / (
                (box_spatial_res ** 2))
    d2T_dz2 = ((sum(z_plus_temps) / len(z_plus_temps)) - (2 * object_temp) + (sum(z_minus_temps) / len(z_minus_temps))) / (
                (box_spatial_res ** 2))
    t_laplacian = (d2T_dx2 + d2T_dy2 + d2T_dz2)  # laplacian is sum of the 2nd derivatives
    dT_dt = matrix_k * t_laplacian
    dT = dT_dt * (delta_time / box_delta_time)  # must normalize delta time to the box delta time
    for vertex_index in vertex_distances:
        if 'C' not in matrix_ids[vertex_index]:
            vertex_distance = vertex_distances[vertex_index]
            distance_ratio = vertex_distance / total_distance
            mesh_temperatures[vertex_index] -= (dT * distance_ratio)
    object_temperatures[object_index] += dT

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
