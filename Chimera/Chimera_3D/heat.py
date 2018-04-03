import numpy as np
import multiprocessing as mp
import time
from . import backends

def conduction(coords, len_coords, x_plus_indices, x_minus_indices, y_plus_indices, y_minus_indices, z_plus_indices,
               z_minus_indices, temperatures, object_ids, spatial_res, conductivities, delta_time, mesh_indices,
               num_workers, multiprocess=False):

    if multiprocess is False:
        update_temps = [0 for i in range(0, len_coords)]  # all of the updated temperatures due to conduction
        dT_dt_list = [0 for i in range(0, len_coords)]
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
        update_temps, dT_dt_list = multiprocess_conduction(coords=coords, len_coords=len_coords, x_plus_indices=x_plus_indices, x_minus_indices=x_minus_indices,
                                        y_plus_indices=y_plus_indices, y_minus_indices=y_minus_indices, z_plus_indices=z_plus_indices,
                                        z_minus_indices=z_minus_indices, temperatures=temperatures, conductivities=conductivities,
                                        spatial_res=spatial_res, delta_time=delta_time,
                                        object_ids=object_ids, mesh_indices=mesh_indices,
                                        num_workers=num_workers, multiprocess=multiprocess)
        return update_temps, dT_dt_list

def multiprocess_conduction(coords, len_coords, x_plus_indices, x_minus_indices, y_plus_indices, y_minus_indices, z_plus_indices,
               z_minus_indices, temperatures, object_ids, spatial_res, conductivities, delta_time, mesh_indices,
               num_workers, multiprocess=True):

    update_temps = []  # all of the updated temperatures due to conduction
    dT_dt_list = []

    jobs = []

    index_chunks = backends.chunk_array(array=mesh_indices, num_chunks=num_workers, len_array=len_coords)
    for process_id, process in enumerate(index_chunks):
        print(process_id, process)
        start_index, stop_index = process[0], process[-1]
        print(start_index, stop_index)
        chunk_coords = coords[start_index:stop_index]
        chunk_x_plus_indices = x_plus_indices[start_index:stop_index]
        chunk_x_minus_indices = x_minus_indices[start_index:stop_index]
        chunk_y_plus_indices = y_plus_indices[start_index:stop_index]
        chunk_y_minus_indices = y_minus_indices[start_index:stop_index]
        chunk_z_plus_indices = z_plus_indices[start_index:stop_index]
        chunk_z_minus_indices = z_minus_indices[start_index:stop_index]
        chunk_temperatures = temperatures[start_index:stop_index]
        chunk_object_ids = object_ids[start_index:stop_index]
        chunk_conductivities = conductivities[start_index:stop_index]
        p = mp.Process(name=str(process_id), target=conduction, args=(chunk_coords, len_coords, chunk_x_plus_indices,
                            chunk_x_minus_indices, chunk_y_plus_indices, chunk_y_minus_indices, chunk_z_plus_indices,
                                                                      chunk_z_minus_indices, chunk_temperatures,
                                                                      chunk_object_ids, spatial_res, chunk_conductivities,
                                                                      delta_time, mesh_indices,
                                                                    num_workers, multiprocess, ))
        jobs.append(p)
        p.start()
    for job in jobs:
        job.join()
            # temp_point = temperatures[index]  # temperature of the specified coordinate position
            # if 'C' not in object_ids[index]:  # makes sure that point z is not a boundary layer
            #     x_plus_index = x_plus_indices[index]  # index of the x+ coordinate position
            #     x_minus_index = x_minus_indices[index]  # index of the x- coordinate position
            #     y_plus_index = y_plus_indices[index]  # index of the y+ coordinate position
            #     y_minus_index = y_minus_indices[index]  # index of the y- coordinate position
            #     z_plus_index = z_plus_indices[index]  # index of the z+ coordinate position
            #     z_minus_index = z_minus_indices[index]  # index of the z- coordinate position
            #     temp_x_plus = temperatures[x_plus_index]  # temperature of the x+ coordinate position
            #     temp_x_minus = temperatures[x_minus_index]  # temperature of the x- coordinate position
            #     temp_y_plus = temperatures[y_plus_index]  # temperature of the y+ coordinate position
            #     temp_y_minus = temperatures[y_minus_index]  # temperature of the y- coordinate position
            #     temp_z_plus = temperatures[z_plus_index]  # temperature of the z+ coordinate position
            #     temp_z_minus = temperatures[z_minus_index]  # temperature of the z- coordinate position
            #     k = conductivities[index]  # conductivities of the material at position of coordinate z
            #
            #     # central difference laplacian is ((f_(x-1) - (2*f_(x)) + f_(x+1)) / (delta_x)^2)
            #     # central difference laplacian is as follows for each vector component
            #     x_temp_laplacian = ((temp_x_plus - (2 * temp_point) + temp_x_minus) / ((spatial_res) ** 2))
            #     y_temp_laplacian = ((temp_y_plus - (2 * temp_point) + temp_y_minus) / ((spatial_res) ** 2))
            #     z_temp_laplacian = ((temp_z_plus - (2 * temp_point) + temp_z_minus) / ((spatial_res) ** 2))
            #     temp_laplacian = x_temp_laplacian + y_temp_laplacian + z_temp_laplacian
            #
            #     # change in temperature with respect to time, dT/dt = -k * laplacian(T)
            #     dT_dt = (k) * temp_laplacian  # the central finite difference heat equation
            #     dT_dt_list.append(dT_dt)
            #
            #     dT = dT_dt * delta_time  # the change in temperature with respect to the finite timestep
            #     new_T = temp_point + dT  # adds dT to the original temperature
            #     update_temps.append(new_T)  # adds the new temperature to the updated temperature list
            # else:  # if it is a boundary layer, it is a fixed temperature
            #     update_temps.append(temp_point)
            #     dT_dt_list.append(0.0)

    return update_temps, dT_dt_list


def viscous_dissipation():
    pass

