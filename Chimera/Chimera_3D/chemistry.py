import pandas as pd
import numpy as np
from statsmodels.formula.api import ols
import os
from pathlib import Path
from . import console, backends


class Chemistry:

    def __init__(self, box):
        self.matrix = {}  # tracks the composition of all components of the matrix
        self.partitioning = {}  # tracks the regression equations of all inserted elements
        self.objects = box.objects  # gets the objects from the box

    def insertMatrixComposition(self, material, composition):
        self.matrix.update({material: composition})
        # automatically calculate the partitioning behavior of the object
        for i in composition:
            if i not in self.partitioning.keys():
                self.regressPartitioning(element=i.lower())
        return None

    def regressPartitioning(self, element):
        data = pd.read_csv(str(Path(__file__).parents[1]) + "/partitioning/{}.csv".format(element.lower()))  # hack for where to get the data for now
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



    def equilibrate(self, object_composition, temperature, pressure, fo2, matrix_material):
        # D = C_solid / C_liquid
        for element in object_composition:
            conc_object = object_composition[element]
            conc_matrix = self.matrix[matrix_material][element]
            print("\n", conc_object, conc_matrix)
            predicted_D = self.partitioning[element]['intercept'] \
                           + (self.partitioning[element]['temperature'] * temperature) \
                           + (self.partitioning[element]['pressure'] * pressure) \
                           + (self.partitioning[element]['fo2'] * fo2)
            current_D = conc_object / conc_matrix
            # we must introduce a factor to change the current partitioning coefficient to the new one
            # adjust = predicted_D / current_D = predicted_D * (conc_liquid / conc_solid)
            adjust = predicted_D / current_D
            if adjust > 1:  # need to increase the concentration in the object
                delta_conc = ((predicted_D * conc_matrix) - conc_object) / (1.0 + predicted_D)
                object_composition[element] += delta_conc
                self.matrix[matrix_material][element] -= delta_conc
            elif adjust < 1:  # need to increase the concentration in the matrix
                delta_conc = ((predicted_D * conc_matrix) - conc_object) / (-1.0 - predicted_D)
                object_composition[element] -= delta_conc
                self.matrix[matrix_material][element] += delta_conc
            else:  # at equilibrium, need to do nothing
                pass
            verify_D = object_composition[element] / self.matrix[matrix_material][element]
            print("Object comp: {}, matrix_comp: {}, adjust: {}, predicted_D: {}, current_D: {}, verify_D: {}".format(
                object_composition[element], self.matrix[matrix_material][element], adjust, predicted_D, current_D,
                verify_D
            ))
        return
