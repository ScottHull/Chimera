import os
import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
from . import console


def temperature_distribution(mesh, verbose=True):

    console.event("Constructing temperature distribution plot...", verbose=verbose)

    t = time.time()
    font = {'size': 10}
    mpl.rc('font', **font)

    objects = mesh['object'].tolist()
    coordinates = mesh['coords'].tolist()
    temperatures = mesh['temperature'].tolist()
    conductivities = mesh['conductivity'].tolist()
    dT_dts = mesh['dT_dt'].tolist()

    fig1 = plt.figure()

    ax1 = fig1.add_subplot(111)
    ax1.plot(coordinates, dT_dts, color='r', linewidth=1.4, linestyle='--')
    ax1.set_xlabel("Depth (m)")
    ax1.set_ylabel("Heat Flux (degK/s)")
    ax1.tick_params('y', colors='r')

    ax2 = ax1.twinx()
    ax2.plot(coordinates, temperatures, color='b', linewidth=2, linestyle='-')
    ax2.set_ylabel("Temperature (degK)")
    ax2.tick_params('y', colors='b')

    fig2 = plt.figure()

    ax3 = fig2.add_subplot(111)
    ax3.plot(coordinates, dT_dts, color='r', linewidth=1.4, linestyle='--')
    ax3.set_xlabel('Depth (m)')
    ax3.set_ylabel('Heat Flux (degK/s)')
    ax3.tick_params('y', colors='r')

    ax4 = ax3.twinx()
    ax4.plot(coordinates, conductivities, color='m', linewidth=2, linestyle='-')
    ax4.set_ylabel('Thermal Conductivity')
    ax4.tick_params('y', colors='m')

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
        color = np.random.rand(3, )
        ax1.axvspan(xmin=min_coord, xmax=max_coord, color=color, alpha=0.2, label=str(object))
        ax3.axvspan(xmin=min_coord, xmax=max_coord, color=color, alpha=0.2, label=str(object))

    ax1.set_title("Temperature Distribution Over Depth")
    ax3.set_title("Thermal Conductivity Over Depth")
    ax1.grid()
    ax1.legend(loc='lower left')
    ax3.grid()
    ax3.legend(loc='lower left')

    console.event("Finished constructing temperature distribution plot! (task took {}s)".format(
        time.time() - t), verbose=verbose)

    plt.show()

def nusselt_distrbution(layers, mesh, verbose=True):

    console.event("Constructing Nusselt distribution plot...", verbose=verbose)

    t = time.time()
    font = {'size': 10}
    mpl.rc('font', **font)


    objects = mesh['object'].tolist()
    coordinates_full = mesh['coords'].tolist()
    dT_dts = mesh['dT_dt'].tolist()
    coordinates = []
    nusselt_nos = []
    coordinates_min_z = layers['min_z'].tolist()
    coordinates_max_z = layers['max_z'].tolist()
    nusselt_nos_list = layers['nusselt'].tolist()
    for index, coord in enumerate(coordinates_min_z):
        coordinates.append(coord)
        coordinates.append(coordinates_max_z[index])
        nusselt_nos.append(nusselt_nos_list[index][0])
        nusselt_nos.append(nusselt_nos_list[index][1])

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(coordinates, nusselt_nos, color='b', linewidth=2, linestyle='-')
    ax1.scatter(coordinates, nusselt_nos, color='b')
    ax1.set_xlabel("Depth (m)")
    ax1.set_ylabel("Nusselt Number")
    ax1.tick_params('y', colors='b')

    ax2 = ax1.twinx()
    ax2.plot(coordinates_full, dT_dts, color='r', linewidth=1.4, linestyle='--')
    ax2.set_ylabel("Heat Flux (degK/s)")
    ax2.tick_params('y', colors='r')

    object_dict = {}
    for index, object in enumerate(objects):
        if object.lower() != 'boundary':
            if object not in object_dict.keys():
                object_dict.update({object: [coordinates_full[index]]})
            else:
                object_dict[object].append(coordinates_full[index])
    for object in object_dict.keys():
        min_coord = min(object_dict[object])
        max_coord = max(object_dict[object])
        color = np.random.rand(3, )
        ax1.axvspan(xmin=min_coord, xmax=max_coord, color=color, alpha=0.2, label=str(object))

    ax1.set_title("Temperature Distribution Over Depth")
    ax1.grid()
    ax1.legend(loc='lower left')

    console.event("Finished constructing Nusselt distribution plot! (task took {}s)".format(
        time.time() - t), verbose=verbose)

    plt.show()
