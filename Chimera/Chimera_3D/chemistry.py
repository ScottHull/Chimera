import pandas as pd
import numpy as np
from statsmodels.formula.api import ols
import os
from pathlib import Path
import copy
from . import console, backends


class Chemistry:

    def __init__(self, box):
        self.box = box
        self.matrix = []  # tracks the composition of all components of the matrix
        self.partitioning = {

        }  # tracks the regression equations of all inserted elements
        self.diffusivities = {

        }
        self.empty_matrix = True

    def insertMatrixComposition(self, index, material, composition, diffusivity):
        if self.empty_matrix:
            self.matrix = [{} for _ in range(len(self.box.mesh['coords']))]
            self.empty_matrix = False
        # automatically calculate the partitioning behavior of the object
        for i in composition:
            self.matrix[index] = composition
            if i.lower() not in self.partitioning.keys():
                self.regressPartitioning(element=i.lower())
            if i.lower() not in self.diffusivities.keys():
                self.diffusivities.update(
                    {
                        i.lower(): diffusivity[i]
                    }
                )
        return None

    def regressPartitioning(self, element):
        data = pd.read_csv(
            str(Path(__file__).parents[1]) + "/partitioning/{}.csv".format(element.lower())
        )  # hack for where to get the data for now
        model = ols("D ~ Temperature + Pressure + fO2", data).fit()
        coeffs = model._results.params
        intercept = coeffs[0]
        temperature_coeff = coeffs[1]
        pressure_coeff = coeffs[2]
        fO2_coeff = coeffs[3]
        self.partitioning.update(
            {
                element:
                    {
                        'intercept': intercept,
                        'temperature': temperature_coeff,
                        'pressure': pressure_coeff,
                        'fo2': fO2_coeff,
                    }
            }
        )
        return self.partitioning

    # def ficks(self, composition, index, coords, x_plus_indices, x_minus_indices, y_plus_indices, y_minus_indices,
    #           z_plus_indices, z_minus_indices, object_ids, len_coords, real_delta_time, update_comps):
    #
    #     # update_comps = [{} for _ in range(0, len_coords)]  # all of the updated compositions due to conduction
    #
    #     delta_time = self.box.delta_time
    #     spatial_res = self.box.spatial_res
    #
    #     for element in composition[index]:
    #         conc_point = composition[index][element]
    #         if 'C' not in object_ids[index]:  # makes sure that point z is not a boundary layer
    #             x_plus_index = x_plus_indices[index]  # index of the x+ coordinate position
    #             x_minus_index = x_minus_indices[index]  # index of the x- coordinate position
    #             y_plus_index = y_plus_indices[index]  # index of the y+ coordinate position
    #             y_minus_index = y_minus_indices[index]  # index of the y- coordinate position
    #             z_plus_index = z_plus_indices[index]  # index of the z+ coordinate position
    #             z_minus_index = z_minus_indices[index]  # index of the z- coordinate position
    #             conc_x_plus = composition[x_plus_index][element]  # composition of the x+ coordinate position
    #             conc_x_minus = composition[x_minus_index][element]  # composition of the x- coordinate position
    #             conc_y_plus = composition[y_plus_index][element]  # composition of the y+ coordinate position
    #             conc_y_minus = composition[y_minus_index][element]  # composition of the y- coordinate position
    #             conc_z_plus = composition[z_plus_index][element]  # composition of the z+ coordinate position
    #             conc_z_minus = composition[z_minus_index][element]  # composition of the z- coordinate position
    #             k = self.diffusivities[element]  # conductivities of the material at position of coordinate z
    #
    #             # central difference laplacian is ((f_(x-1) - (2*f_(x)) + f_(x+1)) / (delta_x)^2)
    #             # central difference laplacian is as follows for each vector component
    #             x_conc_laplacian = ((conc_x_plus - (2 * conc_point) + conc_x_minus) / (spatial_res ** 2))
    #             y_conc_laplacian = ((conc_y_plus - (2 * conc_point) + conc_y_minus) / (spatial_res ** 2))
    #             z_conc_laplacian = ((conc_z_plus - (2 * conc_point) + conc_z_minus) / (spatial_res ** 2))
    #             conc_laplacian = x_conc_laplacian + y_conc_laplacian + z_conc_laplacian
    #
    #             # change in composition with respect to time, dT/dt = -k * laplacian(T)
    #             dC_dt = k * conc_laplacian  # the central finite difference heat equation
    #             # dT_dt_list[index] = dT_dt
    #
    #             dC = dC_dt * (
    #                         delta_time / real_delta_time
    #             )  # the change in composition with respect to the finite normalized timestep
    #             new_C = conc_point + dC  # adds dT to the original composition
    #             update_comps[index][element] = new_C  # adds the new composition to the updated composition list
    #         else:  # if it is a boundary layer, it is a fixed composition
    #             update_comps[index][element] = conc_point
    #             # dT_dt_list[index] = 0.0
    #
    #     return update_comps


    def equilibrate(self, object_concentrations, object_index, vertex_distances, matrix_ids, total_distance,
                    vertex_indices, pressures, temperatures, fO2):

        for element in object_concentrations[object_index]:
            # print(self.matrix)
            # print(object_concentrations)
            cell_matrix_conc = [self.matrix[i][element] for i in vertex_indices]
            cell_pressure = [pressures[i] for i in vertex_indices]
            cell_fO2 = [fO2[i] for i in vertex_indices]
            # D = C_solid / C_liquid
            conc_object = object_concentrations[object_index][element]
            avg_conc_matrix = sum(cell_matrix_conc) / float(len(cell_matrix_conc))
            avg_pressure = sum(cell_pressure) / float(len(cell_pressure))
            avg_fO2 = sum(cell_fO2) / float(len(cell_fO2))
            predicted_D = self.partitioning[element]['intercept'] \
                           + (self.partitioning[element]['temperature'] * temperatures[object_index]) \
                           + (self.partitioning[element]['pressure'] * avg_pressure) \
                           + (self.partitioning[element]['fo2'] * avg_fO2)
            current_D = conc_object / avg_conc_matrix
            adjust = predicted_D / current_D
            if adjust > 1:  # need to increase the concentration in the object
                delta_conc = ((predicted_D * avg_conc_matrix) - conc_object) / (1.0 + predicted_D)
                object_concentrations[object_index][element] += delta_conc
                for vertex_index in vertex_distances:
                    if 'C' not in matrix_ids[vertex_index]:
                        cp_dict = copy.deepcopy(self.matrix[vertex_index])
                        vertex_distance = vertex_distances[vertex_index]
                        distance_ratio = vertex_distance / total_distance
                        cp_dict[element] -= (delta_conc * distance_ratio)
                        self.matrix[vertex_index] = cp_dict
            elif adjust < 1:  # need to increase the concentration in the matrix
                delta_conc = ((predicted_D * avg_conc_matrix) - conc_object) / (-1.0 - predicted_D)
                object_concentrations[object_index][element] -= delta_conc
                for vertex_index in vertex_distances:
                    if 'C' not in matrix_ids[vertex_index]:
                        cp_dict = copy.deepcopy(self.matrix[vertex_index])
                        vertex_distance = vertex_distances[vertex_index]
                        distance_ratio = vertex_distance / total_distance
                        cp_dict[element] += (delta_conc * distance_ratio)
                        self.matrix[vertex_index] = cp_dict
            else:  # at equilibrium, need to do nothing
                pass
        return



    # finite central difference conductivity across entire box if conduction is specified

    # def equilibrate(self, object_composition, temperature, pressure, fo2, matrix_material):
    #     # D = C_solid / C_liquid
    #     for element in object_composition:
    #         conc_object = object_composition[element]
    #         conc_matrix = self.matrix[matrix_material][element]
    #         print("\n", conc_object, conc_matrix)
    #         predicted_D = self.partitioning[element]['intercept'] \
    #                        + (self.partitioning[element]['temperature'] * temperature) \
    #                        + (self.partitioning[element]['pressure'] * pressure) \
    #                        + (self.partitioning[element]['fo2'] * fo2)
    #         current_D = conc_object / conc_matrix
    #         # we must introduce a factor to change the current partitioning coefficient to the new one
    #         # adjust = predicted_D / current_D = predicted_D * (conc_liquid / conc_solid)
    #         adjust = predicted_D / current_D
    #         if adjust > 1:  # need to increase the concentration in the object
    #             delta_conc = ((predicted_D * conc_matrix) - conc_object) / (1.0 + predicted_D)
    #             object_composition[element] += delta_conc
    #             self.matrix[matrix_material][element] -= delta_conc
    #         elif adjust < 1:  # need to increase the concentration in the matrix
    #             delta_conc = ((predicted_D * conc_matrix) - conc_object) / (-1.0 - predicted_D)
    #             object_composition[element] -= delta_conc
    #             self.matrix[matrix_material][element] += delta_conc
    #         else:  # at equilibrium, need to do nothing
    #             pass
    #         verify_D = object_composition[element] / self.matrix[matrix_material][element]
    #         print("Object comp: {}, matrix_comp: {}, adjust: {}, predicted_D: {}, current_D: {}, verify_D: {}".format(
    #             object_composition[element], self.matrix[matrix_material][element], adjust, predicted_D, current_D,
    #             verify_D
    #         ))
    #     return
