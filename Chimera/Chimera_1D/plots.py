import os
import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
from . import console


def temperature_distribution(mesh, verbose=True):
    """
    Plots variables related to the temperature distribution of the linear Chimera model.
    :param mesh:
    :param verbose:
    :return:
    """
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
    """
    Plots variables related to the Nusselt Number distribution of the linear Chimera model.
    :param layers:
    :param mesh:
    :param verbose:
    :return:
    """
    console.event("Constructing Nusselt distribution plot...", verbose=verbose)

    t = time.time()
    font = {'size': 10}
    mpl.rc('font', **font)


    objects = mesh['object'].tolist()
    coordinates_full = mesh['coords'].tolist()
    conductivities = mesh['conductivity'].tolist()
    dT_dts = mesh['dT_dt'].tolist()
    coordinates = []
    nusselt_nos = []
    coordinates_min_z = layers['min_z'].tolist()
    nusselt_nos_list = layers['nusselt'].tolist()
    for index, coord in enumerate(coordinates_min_z):
        coordinates.append(coord)
        nusselt_nos.append(nusselt_nos_list[index])

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(coordinates, nusselt_nos, color='b', linewidth=2, linestyle='-')
    ax1.scatter(coordinates, nusselt_nos, color='b')
    ax1.set_xlabel("Depth (m)")
    ax1.set_ylabel("Nusselt Number")
    ax1.tick_params('y', colors='b')

    ax2 = ax1.twinx()
    ax2.plot(coordinates_full, dT_dts, color='r', linewidth=1.4, linestyle='--')
    ax2.set_ylabel("Heat Flux (degK/s)")
    ax2.tick_params('y', colors='r')

    fig2 = plt.figure()
    ax3 = fig2.add_subplot(111)
    ax3.plot(coordinates, nusselt_nos, color='b', linewidth=2, linestyle='-')
    ax3.scatter(coordinates, nusselt_nos, color='b')
    ax3.set_xlabel("Depth (m)")
    ax3.set_ylabel("Nusselt Number")
    ax3.tick_params('y', colors='b')

    ax4 = ax3.twinx()
    ax4.plot(coordinates_full, conductivities, color='m', linewidth=1.4, linestyle='--')
    ax4.set_ylabel("Thermal Conductivity")
    ax4.tick_params('y', colors='m')

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
        ax3.axvspan(xmin=min_coord, xmax=max_coord, color=color, alpha=0.2, label=str(object))

    ax1.set_title("Temperature Distribution Over Depth")
    ax1.grid()
    ax1.legend(loc='lower left')
    ax3.set_title("Conductivity Distribution Over Depth")
    ax3.grid()
    ax3.legend(loc='lower left')

    console.event("Finished constructing Nusselt distribution plot! (task took {}s)".format(
        time.time() - t), verbose=verbose)

    plt.show()


def nusselt_evolution(df):
    """
    Plot the temporal evolution of the nusselt number for each layer.
    :param df:
    :return:
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for key in df.keys():
        if key != 'time':
            c = np.random.rand(3, )
            ax.plot(df['time'], df[key], color=c, linewidth=1.4, label=key)

    ax.invert_xaxis()
    ax.grid()
    ax.legend(loc='upper right')
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Nusselt Number")
    ax.set_title("Temporal Evolution of the Nusselt Number")

    plt.show()
