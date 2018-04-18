import pandas as pd
import numpy as np
from statsmodels.formula.api import ols
import os
from . import console, backends


class Chemistry:

    def __init__(self, box):
        self.matrix = {}  # tracks the composition of all components of the matrix
        self.partitioning = {}  # tracks the regression equations of all inserted elements
        self.objects = box.objects  # gets the objects from the box

    def insertObjectComposition(self, material, object_id, composition):
        pass

    def insertMatrixComposition(self, material, composition):
        self.matrix.update({material: composition})
        # automatically calculate the partitioning behavior of the object
        for i in composition:
            self.regressPartitioning(element=i)
        return None

    def regressPartitioning(self, element):
        data = pd.read_csv(os.getcwd() + "/{}.csv".format(element.lower()))  # hack for where to get the data for now
        model = ols("partitioncoeff ~ temperature + pressure + fO2", data).fit()
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

    def equilibrate(self, element, temperature, pressure, fo2):
        partitioning = self.partitioning[element]['intercept'] \
                       + self.partitioning[element]['temperature'] * temperature \
                       + self.partitioning[element]['pressure'] * pressure \
                       + self.partitioning[element]['fo2'] * fo2
        return partitioning
