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


    def equilibrate(self, object_concentrations, object_index, vertex_distances, matrix_ids, total_distance,
                    vertex_indices, pressures, temperatures, fO2, spatial_res, object_radius):

        object_volume = (4 / 3) * pi * (object_radius ** 2)
        cell_volume = spatial_res ** 3

        for element in object_concentrations[object_index]:
            cell_matrix_conc = sum([self.matrix[i][element] for i in vertex_indices]) / cell_volume
            cell_pressure = [pressures[i] for i in vertex_indices]
            cell_fO2 = [fO2[i] for i in vertex_indices]
            # D = C_solid / C_liquid
            conc_object = object_concentrations[object_index][element] / object_volume
            avg_pressure = sum(cell_pressure) / float(len(cell_pressure))
            avg_fO2 = sum(cell_fO2) / float(len(cell_fO2))
            predicted_D = self.partitioning[element]['intercept'] \
                          + (self.partitioning[element]['temperature'] * temperatures[object_index]) \
                          + (self.partitioning[element]['pressure'] * avg_pressure) \
                          + (self.partitioning[element]['fo2'] * avg_fO2)
            current_D = conc_object / cell_matrix_conc
            adjust = predicted_D / current_D

            if adjust > 1:  # need to increase the concentration in the object
                delta_conc = ((predicted_D * cell_matrix_conc) - conc_object) / (1.0 + predicted_D)
                object_concentrations[object_index][element] += (delta_conc * object_volume)
                for vertex_index in vertex_indices:
                    if 'C' not in matrix_ids[vertex_index]:
                        cp_dict = copy.deepcopy(self.matrix[vertex_index])
                        cp_dict[element] -= \
                            delta_conc * (cp_dict[element] / (cell_matrix_conc * cell_volume)) * cell_volume
                        self.matrix[vertex_index] = cp_dict

            elif adjust < 1:  # need to increase the concentration in the matrix
                delta_conc = ((predicted_D * cell_matrix_conc) - conc_object) / (-1.0 - predicted_D)
                object_concentrations[object_index][element] -= (delta_conc * object_volume)
                for vertex_index in vertex_indices:
                    if 'C' not in matrix_ids[vertex_index]:
                        cp_dict = copy.deepcopy(self.matrix[vertex_index])
                        cp_dict[element] += \
                            delta_conc * (cp_dict[element] / (cell_matrix_conc * cell_volume)) * cell_volume
                        self.matrix[vertex_index] = cp_dict

            else:  # at equilibrium, need to do nothing
                pass

            # print("\nADJUST: {}, OBJECT_CONC: {}, LIQUID_CONC: {}, PREDICTED_D: {}, CONFIRM_D: {}".format(
            #     adjust, object_concentrations[object_index][element],
            #     sum([self.matrix[i][element] for i in vertex_indices]),
            #     predicted_D,
            #     (object_concentrations[object_index][element] / object_volume) /
            #     (sum([self.matrix[i][element] for i in vertex_indices]) / cell_volume)
            # ))

        return

    def resetMatrixComp(self, new_matrix_comp):
        self.matrix = new_matrix_comp
        return