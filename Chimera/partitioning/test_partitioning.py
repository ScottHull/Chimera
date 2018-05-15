import pandas as pd
import numpy as np
from statsmodels.formula.api import ols
import os
import matplotlib.pyplot as plt

def regressPartitioning(element):
    partitioning = {}
    data = pd.read_csv(os.getcwd() + "/{}.csv".format(
        element.lower()))  # hack for where to get the data for now
    model = ols("D ~ Temperature + Pressure + fO2", data).fit()
    coeffs = model._results.params
    intercept = coeffs[0]
    temperature_coeff = coeffs[1]
    pressure_coeff = coeffs[2]
    fO2_coeff = coeffs[3]
    partitioning.update({element:
        {
            'intercept': intercept,
            'temperature': temperature_coeff,
            'pressure': pressure_coeff,
            'fo2': fO2_coeff,
        }
    })
    return partitioning

def regressPartitioning_fO2_D(element):
    data = pd.read_csv(os.getcwd() + "/{}.csv".format(
        element.lower()))  # hack for where to get the data for now
    model = ols("D ~ fO2", data).fit()
    coeffs = model._results.params
    intercept = coeffs[0]
    fO2_coeff = coeffs[1]
    return fO2_coeff, intercept

element = 'w'
fO2_coeff, intercept = regressPartitioning_fO2_D(element=element)
data = pd.read_csv(os.getcwd() + "/{}.csv".format(
        element.lower()))

fig = plt.figure()
fO2_D = fig.add_subplot(111)
fO2_D.plot(data['fO2'], [(fO2_coeff * i) + intercept for i in data['fO2']])
fO2_D.scatter(data['fO2'], data['D'])
plt.show()

