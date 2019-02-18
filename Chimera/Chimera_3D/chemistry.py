import pandas as pd
import numpy as np
from statsmodels.formula.api import ols
import os
from pathlib import Path
import copy
from math import pi
from . import console, backends


class Chemistry:

    def __init__(self, box):
        self.box = box
        self.matrix = []  # tracks the composition of all components of the matrix
        self.partitioning = {}  # tracks the regression equations of all inserted elements
        self.diffusivities = {}
        self.models = {}
        self.empty_matrix = True
        self.track_distribution_coeffs = {
            "z-depth": [],
            'object': [],
            'cell_temperature': [],
            'cell_pressure': [],
            'cell_fO2': [],
            'cell_conc': [],
            'object_conc': [],
            'object_temperature': [],
            'D': []
        }

    def insertMatrixComposition(self, index, material, composition, diffusivity):

        if self.empty_matrix:
            self.matrix = [{} for _ in range(len(self.box.mesh['coords']))]
            self.empty_matrix = False
        # automatically calculate the partitioning behavior of the object
        for i in composition:
            self.matrix[index] = composition
            if i.lower() not in self.diffusivities.keys():
                self.diffusivities.update(
                    {
                        i.lower(): diffusivity[i]
                    }
                )
        return None


    # def insertModel(self, element, alpha, beta, chi, delta, epsilon):
    #
    #     # eg. Cottrell et al. (2009) metal-silicate W partitioning model:
    #     # log(D) = alpha + beta * (delta IW) + chi * (nbo/t) + delta * (1/T) + epsilon(P/T)
    #
    #     self.partitioning.update({
    #         element: {
    #             "alpha": alpha,
    #             "beta": beta,
    #             "chi": chi,
    #             "delta": delta,
    #             "epsilon": epsilon,
    #         }
    #     })
    #
    #     return None


    def applyModel(self, element, pressure, temperature, fO2, nbo_t=2.6):

        # eg. Cottrell et al. (2009) metal-silicate W partitioning model:
        # log(D) = alpha + beta * (delta IW) + chi * (nbo/t) + delta * (1/T) + epsilon(P/T)

        # alpha = self.partitioning[element]['alpha']
        # beta = self.partitioning[element]['beta']
        # chi = self.partitioning[element]['chi']
        # delta = self.partitioning[element]['delta']
        # epsilon = self.partitioning[element]['epsilon']

        coeffs = {
            'alpha': 0,
            'beta': 0,
            'chi': 0,
            'delta': 0,
            'epsilon': 0
        }

        if 0.0 <= pressure <= 2:
            coeffs['alpha'] = 1.11
            coeffs['beta'] = -1.18
            coeffs['chi'] = -0.85
            coeffs['delta'] = 1680
            coeffs['epsilon'] = 487
        # elif pressure == 2:
        #     if 2300 < temperature < 2600:
        #         coeffs['alpha'] = 0.84
        #         coeffs['beta'] = -1.22
        #         coeffs['chi'] = -0.85
        #         coeffs['delta'] = 3245
        #         coeffs['epsilon'] = 487
        # elif pressure == 6:
        #     if 2300 < temperature < 2700:
        #         coeffs['alpha'] = 1.17
        #         coeffs['beta'] = -1.06
        #         coeffs['chi'] = -0.90
        #         coeffs['delta'] = 3337
        #         coeffs['epsilon'] = 487
        # elif 2 < pressure:
        #     coeffs['alpha'] = 1.05
        #     coeffs['beta'] = -1.10
        #     coeffs['chi'] = -0.84
        #     coeffs['delta'] = 3588
        #     coeffs['epsilon'] = -102
        #
        # if 0.5 <= pressure <= 2:
        #     if 2100 < temperature < 2600:
        #         coeffs['alpha'] = 1.11
        #         coeffs['beta'] = -1.18
        #         coeffs['chi'] = -0.85
        #         coeffs['delta'] = 1680
        #         coeffs['epsilon'] = 487
        # elif pressure == 2:
        #     if 2300 < temperature < 2600:
        #         coeffs['alpha'] = 0.84
        #         coeffs['beta'] = -1.22
        #         coeffs['chi'] = -0.85
        #         coeffs['delta'] = 3245
        #         coeffs['epsilon'] = 487
        # elif pressure == 6:
        #     if 2300 < temperature < 2700:
        #         coeffs['alpha'] = 1.17
        #         coeffs['beta'] = -1.06
        #         coeffs['chi'] = -0.90
        #         coeffs['delta'] = 3337
        #         coeffs['epsilon'] = 487
        # elif 2 < pressure < 18:
        #     if 2300 < temperature < 2700:
        #         coeffs['alpha'] = 1.05
        #         coeffs['beta'] = -1.10
        #         coeffs['chi'] = -0.84
        #         coeffs['delta'] = 3588
        #         coeffs['epsilon'] = -102

        alpha = coeffs['alpha']
        beta = coeffs['beta']
        chi = coeffs['chi']
        delta = coeffs['delta']
        epsilon = coeffs['epsilon']

        logD = alpha + (beta * fO2) + (chi * nbo_t) + (delta * (1/temperature)) + (epsilon * (pressure/temperature))
        D = 10**logD

        return D


    # def regressPartitioning(self, element):
    #
    #     # hack for where to get the data for now
    #     data = pd.read_csv(str(Path(__file__).parents[1]) + "/partitioning/{}.csv".format(element.lower()))
    #     model = ols("D ~ Temperature + Pressure + fO2", data=data).fit()
    #
    #     self.models.update(
    #         {
    #             element: {
    #                 model.summary()
    #             }
    #         }
    #     )
    #
    #     coeffs = model._results.params
    #     intercept = coeffs[0]
    #     temperature_coeff = coeffs[1]
    #     pressure_coeff = coeffs[2]
    #     fO2_coeff = coeffs[3]
    #
    #     self.partitioning.update(
    #         {
    #             element: {
    #                 'intercept': intercept,
    #                 'temperature': temperature_coeff,
    #                 'pressure': pressure_coeff,
    #                 'fo2': fO2_coeff,
    #             }
    #         }
    #     )
    #     return self.partitioning


    def equilibrate(self, object_moles, object_index, vertex_distances, matrix_ids, total_distance,
                     vertex_indices, pressures, temperatures, object_temperatures, fO2, spatial_res, object_radius,
                    z_depth, object_id):

        object_volume = (4 / 3) * pi * (object_radius ** 3)
        cell_volume = spatial_res ** 3

        for element in object_moles[object_index]:
            moles_cell = sum([self.matrix[i][element] for i in vertex_indices])
            cell_matrix_conc = moles_cell / cell_volume
            cell_pressure = [pressures[i] for i in vertex_indices]
            cell_temperature = [temperatures[i] for i in vertex_indices]
            cell_fO2 = [fO2[i] for i in vertex_indices]
            moles_object = object_moles[object_index][element]
            conc_object = moles_object / object_volume
            avg_pressure = sum(cell_pressure) / float(len(cell_pressure))
            avg_fO2 = sum(cell_fO2) / float(len(cell_fO2))
            avg_temperature = sum(cell_temperature) / float(len(cell_temperature))
            object_temperature = object_temperatures[object_index]
            predicted_D = self.applyModel(
                element=element,
                temperature=avg_temperature,
                pressure=avg_pressure,
                fO2=avg_fO2,
                nbo_t=2.6
            )
            current_D = conc_object / cell_matrix_conc
            # print("CURRENT D ",current_D)
            adjust = predicted_D / current_D

            if adjust > 1.0:
                adj_matrix = (moles_cell /
                                (1 + (3 * ((spatial_res)**3) * ((4 * pi * (object_radius**3) * predicted_D)**(-1)))))
                adj_object = (moles_object /
                              (1 + (4 * pi * (object_radius**3) * predicted_D) * ((3**(-1)) * (spatial_res**(-3)))))
                adj_moles = adj_matrix - adj_object

                # adjust the moles of the element in the object and matrix, respectively
                object_moles[object_index][element] += adj_moles
                for i in vertex_indices:
                    copy_matrix_dict = copy.deepcopy(self.matrix[i])
                    copy_matrix_dict[element] -= (adj_moles / len(vertex_indices))
                    self.matrix[i] = copy_matrix_dict

                new_conc_object = object_moles[object_index][element] / object_volume
                new_conc_cell = sum([self.matrix[i][element] for i in vertex_indices]) / cell_volume
                new_D = new_conc_object / new_conc_cell

                self.track_distribution_coeffs['object_temperature'].append(object_temperature)
                self.track_distribution_coeffs['cell_conc'].append(new_conc_object)

            elif adjust < 1.0:
                adj_matrix = (moles_cell /
                                ((3 * ((spatial_res)**3) * ((4 * pi * (object_radius**3) * predicted_D)**(-1))) - 1))
                adj_object = (moles_object /
                              (1 - (4 * pi * (object_radius ** 3) * predicted_D) * ((3 ** (-1)) * (spatial_res ** (-3)))))
                adj_moles = adj_object - adj_matrix

                # adjust the moles of the element in the object and matrix, respectively
                object_moles[object_index][element] -= adj_moles
                for i in vertex_indices:
                    copy_matrix_dict = copy.deepcopy(self.matrix[i])
                    copy_matrix_dict[element] += (adj_moles / len(vertex_indices))
                    self.matrix[i] = copy_matrix_dict

                new_conc_object = object_moles[object_index][element] / object_volume
                new_conc_cell = sum([self.matrix[i][element] for i in vertex_indices]) / cell_volume
                new_D = new_conc_object / new_conc_cell

                self.track_distribution_coeffs['object_temperature'].append(object_temperature)
                self.track_distribution_coeffs['cell_conc'].append(new_conc_object)

            else:
                pass

            self.track_distribution_coeffs['z-depth'].append(z_depth)
            self.track_distribution_coeffs['object'].append(object_id)
            self.track_distribution_coeffs['cell_temperature'].append(avg_temperature)
            self.track_distribution_coeffs['cell_pressure'].append(avg_pressure)
            self.track_distribution_coeffs['cell_fO2'].append(avg_fO2)
            self.track_distribution_coeffs['cell_conc'].append(new_conc_cell)
            self.track_distribution_coeffs['D'].append(predicted_D)




    # def equilibrate(self, object_concentrations, object_index, vertex_distances, matrix_ids, total_distance,
    #                 vertex_indices, pressures, temperatures, fO2, spatial_res, object_radius):
    #
    #     object_volume = (4 / 3) * pi * (object_radius ** 2)
    #     cell_volume = spatial_res ** 3
    #
    #     for element in object_concentrations[object_index]:
    #         cell_matrix_conc = sum([self.matrix[i][element] for i in vertex_indices]) / cell_volume
    #         cell_pressure = [pressures[i] for i in vertex_indices]
    #         cell_fO2 = [fO2[i] for i in vertex_indices]
    #         # D = C_solid / C_liquid
    #         conc_object = object_concentrations[object_index][element] / object_volume
    #         avg_pressure = sum(cell_pressure) / float(len(cell_pressure))
    #         avg_fO2 = sum(cell_fO2) / float(len(cell_fO2))
    #         predicted_D = self.partitioning[element]['intercept'] \
    #                       + (self.partitioning[element]['temperature'] * temperatures[object_index]) \
    #                       + (self.partitioning[element]['pressure'] * avg_pressure) \
    #                       + (self.partitioning[element]['fo2'] * avg_fO2)
    #         current_D = conc_object / cell_matrix_conc
    #         adjust = predicted_D / current_D
    #
    #         if adjust > 1:  # need to increase the concentration in the object
    #             delta_conc = ((predicted_D * cell_matrix_conc) - conc_object) / (1.0 + predicted_D)
    #             object_concentrations[object_index][element] += (delta_conc * object_volume)
    #             for vertex_index in vertex_indices:
    #                 if 'C' not in matrix_ids[vertex_index]:
    #                     cp_dict = copy.deepcopy(self.matrix[vertex_index])
    #                     cp_dict[element] -= \
    #                         delta_conc * (cp_dict[element] / (cell_matrix_conc * cell_volume)) * cell_volume
    #                     self.matrix[vertex_index] = cp_dict
    #
    #         elif adjust < 1:  # need to increase the concentration in the matrix
    #             delta_conc = ((predicted_D * cell_matrix_conc) - conc_object) / (-1.0 - predicted_D)
    #             object_concentrations[object_index][element] -= (delta_conc * object_volume)
    #             for vertex_index in vertex_indices:
    #                 if 'C' not in matrix_ids[vertex_index]:
    #                     cp_dict = copy.deepcopy(self.matrix[vertex_index])
    #                     cp_dict[element] += \
    #                         delta_conc * (cp_dict[element] / (cell_matrix_conc * cell_volume)) * cell_volume
    #                     self.matrix[vertex_index] = cp_dict
    #
    #         else:  # at equilibrium, need to do nothing
    #             pass
    #
    #         # print("\nADJUST: {}, OBJECT_CONC: {}, LIQUID_CONC: {}, PREDICTED_D: {}, CONFIRM_D: {}".format(
    #         #     adjust, object_concentrations[object_index][element],
    #         #     sum([self.matrix[i][element] for i in vertex_indices]),
    #         #     predicted_D,
    #         #     (object_concentrations[object_index][element] / object_volume) /
    #         #     (sum([self.matrix[i][element] for i in vertex_indices]) / cell_volume)
    #         # ))
    #
    #     return

    def resetMatrixComp(self, new_matrix_comp):

        self.matrix = new_matrix_comp
        return