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
                self.regressPartitioning(element=i)
        return None

    def regressPartitioning(self, element):
        data = pd.read_csv(str(Path(__file__).parents[1]) + "/partitioning/{}.csv".format(element.lower()))  # hack for where to get the data for now
        model = ols("D ~ Temperature + Pressure + fO2", data).fit()
        coeffs = model._results.params
        intercept = coeffs[0]
        temperature_coeff = coeffs[1]
        pressure_coeff = coeffs[2]
        fO2_coeff = coeffs[3]
        self.partitioning.update({element:
                                      {
                                          'intercept': intercept,
                                            'temperature': temperature_coeff,
                                            'pressure': pressure_coeff,
                                            'fo2': fO2_coeff,
                                            }
        })
        return self.partitioning



    def equilibrate(self, object_composition, temperature, pressure, fo2, matrix_composition):
        # D = C_liquid / C_solid
        for element in object_composition:
            D = self.partitioning[element]['intercept'] \
                           + self.partitioning[element]['temperature'] * temperature \
                           + self.partitioning[element]['pressure'] * pressure \
                           + self.partitioning[element]['fo2'] * fo2
            new_object_conc = object_composition[element] + (matrix_composition[element] / D)  # c_solid = C_liquid / D
            new_matrix_conc = object_composition[element] * D # C_liquid = D * C_solid
            return new_object_conc, new_matrix_conc
