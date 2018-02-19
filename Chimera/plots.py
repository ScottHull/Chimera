import os
import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
from . import console



def temperature_distribution(df):
    ax = plt.figure().add_subplot(111)

    objects = df['object'].tolist()
    coordinates = df['coords'].tolist()
    temperatures = df['temperature'].tolist()

    ax.plot(coordinates, temperatures, linewidth=2)
    ax.set_xlabel("Depth (m)")
    ax.set_ylabel("Temperature (degK)")
    ax.set_title("Temperature Distribution Over Depth")
    ax.grid()

    object_dict = {}
    for index, object in enumerate(objects):
        if object.lower() != 'boundary':
            if object not in object_dict.keys():
                object_dict.update({object: [coordinates[index]]})
            else:
                object_dict[object].append(coordinates[index])

    for object in object_dict.keys():
        min_coord = min(object_dict[object])
        max_coord = max(object_dict[object])
        ax.axvspan(xmin=min_coord, xmax=max_coord, color=np.random.rand(3,), alpha=0.5, label=str(object))

    ax.legend(loc='upper right')

    plt.show()

