import os
import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
from . import console


def plot_all(mesh_df, layers_df, nusselt_df, verbose=True, save=False, show=True):
    temperature_distribution(mesh_df=mesh_df, verbose=verbose, save=save, show=show)
    nusselt_distrbution(layers_df=layers_df, mesh_df=mesh_df, verbose=verbose, save=save, show=show)
    nusselt_evolution(nusselt_df=nusselt_df, verbose=verbose, save=save, show=show)

def temperature_distribution(mesh_df, verbose=True, save=False, show=True):
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

    objects = mesh_df['object'].tolist()
    coordinates = mesh_df['coords'].tolist()
    temperatures = mesh_df['temperature'].tolist()
    conductivities = mesh_df['conductivity'].tolist()
    dT_dts = mesh_df['dT_dt'].tolist()

    fig1 = plt.figure(figsize=(8.0, 5.0))  # depth vs heat flux, depth vs temperature

    ax1 = fig1.add_subplot(111)
    ax1.plot(coordinates, dT_dts, color='r', linewidth=1.4, linestyle='--')
    ax1.set_xlabel("Depth (m)")
    ax1.set_ylabel("Heat Flux (degK/s)")
    ax1.tick_params('y', colors='r')

    ax2 = ax1.twinx()
    ax2.plot(coordinates, temperatures, color='b', linewidth=2, linestyle='-')
    ax2.set_ylabel("Temperature (degK)")
    ax2.tick_params('y', colors='b')

    fig2 = plt.figure(figsize=(8.0, 5.0))  # depth vs heat flux, depth vs thermal conductivity

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

    if show is True:
        plt.show()
    if save is True:
        fig1.tight_layout()
        fig2.tight_layout()
        fig1_name = "temp_distrib_fig1.png"
        fig2_name = "temp_distrib_fig2.png"
        if fig1_name in os.listdir(os.getcwd()):
            os.remove(fig1_name)
        if fig2_name in os.listdir(os.getcwd()):
            os.remove(fig2_name)
        fig1.savefig(fig1_name, format='png')
        fig2.savefig(fig2_name, format='png')

def nusselt_distrbution(layers_df, mesh_df, verbose=True, save=False, show=True):
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


    objects = mesh_df['object'].tolist()
    coordinates_full = mesh_df['coords'].tolist()
    conductivities = mesh_df['conductivity'].tolist()
    dT_dts = mesh_df['dT_dt'].tolist()
    coordinates = []
    nusselt_nos = []
    coordinates_min_z = layers_df['min_z'].tolist()
    nusselt_nos_list = layers_df['nusselt'].tolist()
    for index, coord in enumerate(coordinates_min_z):
        coordinates.append(coord)
        nusselt_nos.append(nusselt_nos_list[index])

    fig1 = plt.figure(figsize=(8.0, 5.0))  # depth vs nusselt number, depth vs heat flux

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

    fig2 = plt.figure(figsize=(8.0, 5.0))  # depth vs nusselt number, depth vs thermal conductivity

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

    ax1.set_title("Nusselt No. Distribution Over Depth")
    ax1.grid()
    ax1.legend(loc='lower left')
    ax3.set_title("Nusselt No. Distribution Over Depth")
    ax3.grid()
    ax3.legend(loc='lower left')

    console.event("Finished constructing Nusselt distribution plot! (task took {}s)".format(
        time.time() - t), verbose=verbose)

    if show is True:
        plt.show()
    if save is True:
        fig1.tight_layout()
        fig2.tight_layout()
        fig1_name = "nusselt_distrib_fig1.png"
        fig2_name = "nusselt_distrib_fig2.png"
        if fig1_name in os.listdir(os.getcwd()):
            os.remove(fig1_name)
        if fig2_name in os.listdir(os.getcwd()):
            os.remove(fig2_name)
        fig1.savefig(fig1_name, format='png')
        fig2.savefig(fig2_name, format='png')


def nusselt_evolution(nusselt_df, verbose=True, save=False, show=True):
    """
    Plot the temporal evolution of the nusselt number for each layer.
    :param df:
    :return:
    """
    fig = plt.figure(figsize=(8.0, 5.0))  # time vs nusselt number evolution

    ax = fig.add_subplot(111)

    for key in nusselt_df.keys():
        if key != 'time':
            c = np.random.rand(3, )
            ax.plot(nusselt_df['time'], nusselt_df[key], color=c, linewidth=1.4, label=key)

    ax.invert_xaxis()
    ax.grid()
    ax.legend(loc='upper right')
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Nusselt Number")
    ax.set_title("Temporal Evolution of the Nusselt Number")

    if show is True:
        plt.show()
    if save is True:
        fig.tight_layout()
        fig_name = "nusselt_evolution.png"
        if fig_name in os.listdir(os.getcwd()):
            os.remove(fig_name)
        fig.savefig(fig_name, format='png')
