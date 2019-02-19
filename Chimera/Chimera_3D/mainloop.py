import numpy as np
import copy
from .backends import override_timestep

def modelLoop(conduction, chemical_diffusion, coords, chemistry, len_coords, x_plus_indices, x_minus_indices, y_plus_indices,
              y_minus_indices, z_plus_indices, z_minus_indices, temperatures, object_ids, spatial_res, spatial_sigfigs,
              conductivities, delta_time, therm_diffusivities, mesh_indices, num_workers, verbose,
              multiprocess=False):

    # the thermal & chemical timestep relative to the possible user defined timestep
    therm_real_delta_time = override_timestep(
            timestep=False,
            conductivities=list(conductivities),
            spatial_res=spatial_res,
            spatial_sigfigs=spatial_sigfigs,
            diffusivities=therm_diffusivities,
            verbose=verbose
        )
    chem_real_delta_time = override_timestep(
            timestep=False,
            conductivities=list(conductivities),
            spatial_res=spatial_res,
            spatial_sigfigs=spatial_sigfigs,
            diffusivities=np.array([chemistry.diffusivities[i] for i in chemistry.diffusivities]),
            verbose=verbose
        )
    update_temps = [0.0 for _ in
                    range(0, len_coords)]  # all of the updated temperatures due to conduction
    update_dT_dt = [0.0 for _ in range(0, len_coords)]

    update_comps = copy.deepcopy(chemistry.matrix)  # all of the updated compositions due to conduction

    for index, coord in enumerate(coords):
        x_plus_index = x_plus_indices[index]  # index of the x+ coordinate position
        x_minus_index = x_minus_indices[index]  # index of the x- coordinate position
        y_plus_index = y_plus_indices[index]  # index of the y+ coordinate position
        y_minus_index = y_minus_indices[index]  # index of the y- coordinate position
        z_plus_index = z_plus_indices[index]  # index of the z+ coordinate position
        z_minus_index = z_minus_indices[index]  # index of the z- coordinate position

        # fourier's law of thermal conduction
        if conduction:
            temp_point = temperatures[index]  # temperature of the specified coordinate position
            if 'C' not in object_ids[index]:  # makes sure that point z is not a boundary layer
                temp_x_plus = temperatures[x_plus_index]  # temperature of the x+ coordinate position
                temp_x_minus = temperatures[x_minus_index]  # temperature of the x- coordinate position
                temp_y_plus = temperatures[y_plus_index]  # temperature of the y+ coordinate position
                temp_y_minus = temperatures[y_minus_index]  # temperature of the y- coordinate position
                temp_z_plus = temperatures[z_plus_index]  # temperature of the z+ coordinate position
                temp_z_minus = temperatures[z_minus_index]  # temperature of the z- coordinate position
                k = therm_diffusivities[index]  # conductivities of the material at position of coordinate z

                # central difference laplacian is ((f_(x-1) - (2*f_(x)) + f_(x+1)) / (delta_x)^2)
                # central difference laplacian is as follows for each vector component
                x_temp_laplacian = ((temp_x_plus - (2 * temp_point) + temp_x_minus) / ((spatial_res) ** 2))
                y_temp_laplacian = ((temp_y_plus - (2 * temp_point) + temp_y_minus) / ((spatial_res) ** 2))
                z_temp_laplacian = ((temp_z_plus - (2 * temp_point) + temp_z_minus) / ((spatial_res) ** 2))
                temp_laplacian = x_temp_laplacian + y_temp_laplacian + z_temp_laplacian

                # change in temperature with respect to time, dT/dt = -k * laplacian(T)
                dT_dt = k * temp_laplacian  # the central finite difference heat equation
                update_dT_dt[index] = dT_dt

                dT = dT_dt * (
                            delta_time / therm_real_delta_time)  # the change in temperature with respect to the finite normalized timestep
                new_T = temp_point + dT  # adds dT to the original temperature
                update_temps[index] = new_T  # adds the new temperature to the updated temperature list
            else:  # if it is a boundary layer, it is a fixed temperature
                update_temps[index] = temp_point
                update_dT_dt[index] = 0.0

        # fick's 2nd law of chemical diffusion
        if chemical_diffusion:
            for element in chemistry.matrix[index]:
                conc_point = chemistry.matrix[index][element]
                if 'C' not in object_ids[index]:  # makes sure that point z is not a boundary layer
                    conc_x_plus = chemistry.matrix[x_plus_index][element]  # concentration of the x+ coordinate position
                    conc_x_minus = chemistry.matrix[x_minus_index][element]  # concentration of the x- coordinate position
                    conc_y_plus = chemistry.matrix[y_plus_index][element]  # concentration of the y+ coordinate position
                    conc_y_minus = chemistry.matrix[y_minus_index][element]  # concentration of the y- coordinate position
                    conc_z_plus = chemistry.matrix[z_plus_index][element]  # concentration of the z+ coordinate position
                    conc_z_minus = chemistry.matrix[z_minus_index][element]  # concentration of the z- coordinate position
                    k = chemistry.diffusivities[element]  # chemical diffusivities of the material at position of coordinate z

                    # central difference laplacian is ((f_(x-1) - (2*f_(x)) + f_(x+1)) / (delta_x)^2)
                    # central difference laplacian is as follows for each vector component
                    x_conc_laplacian = ((conc_x_plus - (2 * conc_point) + conc_x_minus) / (spatial_res ** 2))
                    y_conc_laplacian = ((conc_y_plus - (2 * conc_point) + conc_y_minus) / (spatial_res ** 2))
                    z_conc_laplacian = ((conc_z_plus - (2 * conc_point) + conc_z_minus) / (spatial_res ** 2))
                    conc_laplacian = x_conc_laplacian + y_conc_laplacian + z_conc_laplacian

                    # change in chemistry.matrix with respect to time, dT/dt = -k * laplacian(T)
                    dC_dt = k * conc_laplacian  # the central finite difference heat equation
                    # dT_dt_list[index] = dT_dt

                    # the change in chemistry.matrix with respect to the finite normalized timestep
                    dC = dC_dt * (
                            delta_time / chem_real_delta_time
                    )
                    # adds dT to the original chemistry.matrix
                    new_C = conc_point + dC
                    # adds the new chemistry.matrix to the updated chemistry.matrix list
                    update_comps[index][element] = new_C
                # if it is a boundary layer, it is a fixed chemistry.matrix
                else:
                    update_comps[index][element] = conc_point


    return update_temps, update_dT_dt, update_comps

def noBoundaryChemicalDiffusion(
        conduction, chem, coords, chemistry, len_coords, x_plus_indices, x_minus_indices, y_plus_indices,
        y_minus_indices, z_plus_indices, z_minus_indices, temperatures, object_ids, spatial_res,
        spatial_sigfigs,
        conductivities, delta_time, therm_diffusivities, mesh_indices, num_workers, verbose,
        multiprocess=False):

        # the thermal & chemical timestep relative to the possible user defined timestep
        therm_real_delta_time = override_timestep(
            timestep=False,
            conductivities=list(conductivities),
            spatial_res=spatial_res,
            spatial_sigfigs=spatial_sigfigs,
            diffusivities=therm_diffusivities,
            verbose=verbose
        )
        chem_real_delta_time = override_timestep(
            timestep=False,
            conductivities=list(conductivities),
            spatial_res=spatial_res,
            spatial_sigfigs=spatial_sigfigs,
            diffusivities=np.array([chemistry.diffusivities[i] for i in chemistry.diffusivities]),
            verbose=verbose
        )
        update_temps = [0.0 for _ in
                        range(0, len_coords)]  # all of the updated temperatures due to conduction
        update_dT_dt = [0.0 for _ in range(0, len_coords)]

        update_comps = copy.deepcopy(chemistry.matrix)  # all of the updated compositions due to conduction

        for index, coord in enumerate(coords):
            x_plus_index = x_plus_indices[index]  # index of the x+ coordinate position
            x_minus_index = x_minus_indices[index]  # index of the x- coordinate position
            y_plus_index = y_plus_indices[index]  # index of the y+ coordinate position
            y_minus_index = y_minus_indices[index]  # index of the y- coordinate position
            z_plus_index = z_plus_indices[index]  # index of the z+ coordinate position
            z_minus_index = z_minus_indices[index]  # index of the z- coordinate position
            object_id_z_plus = object_ids[z_plus_index]
            object_id_z_minus = object_ids[z_minus_index]

            # check if the model loop is near a boundary layer
            if "C" in object_id_z_minus:
                pass
            elif "C" in object_id_z_plus:
                pass
