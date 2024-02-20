# Command to launch the test :
# python3 display_curves.py ../../../build/tests/geometryRTheta/advection_2d_rp/<selected advection test>

"""
Plot the advected function at regular different time steps.

Parameters
----------
executable : string
    Path to the executable of the advection simulation.

rmin : float
    Minimum value of the r values.
    Possibiliy for the user to change it by addind the two
    new sizes in the command line.

rmax : float
    Maximum value of the r values.
    Possibiliy for the user to change it by addind the two
    new sizes in the command line.

Nr : int
    Number of break points in r dimension
    Possibiliy for the user to change it by addind the two
    new sizes in the command line.

Nth : int
    Number of break points in theta dimension
    Possibiliy for the user to change it by addind the two
    new sizes in the command line.

dt : float
    Time step.
    Possibiliy for the user to change it by addind the two
    new sizes in the command line.

T : float
    Final time.
    Possibiliy for the user to change it by addind the two
    new sizes in the command line.

curves : bool
    Boolean to select if the values of the advected function are saved.
    Possibiliy for the user to change it by addind the two
    new sizes in the command line.

feet : bool
    Boolean to select if the characteristic feet are saved.
    Possibiliy for the user to change it by addind the two
    new sizes in the command line.
"""

import os
import numpy as np
import matplotlib.pyplot as plt


from advection_functions import set_input, execute, treatment, get_simulation_config

# Definition of functions ---------------------------------------
def set_axis(ax, x_label, y_label, z_label, ax_label):
    """
    Set the labels of the axis of the plot.

    Parameters
    ----------
    ax: matplotlib.axis
        The axis object.
    x_label: string
        The label for the x axis.
    y_label: string
        The label for the y axis.
    z_label: string
        The label for the z axis.
    ax_label: string
        To set property on the axis.
    """
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    ax.axis(ax_label)


def set_surface(fig, nb_x, nb_y, n_idx, X, Y, Z, elev, azim, roll):
    """
    Set a plot_trisurf curve.

    Parameters
    ----------
    fig: matplotlib.figure
        The figure object of the plot.
    nb_x: int
        The number of rows of the figure.
    nb_y: int
        The number of columns of the figure.
    n_idx: int
        The number of the selected subplot.
    X: array
        The coordinates in the x dimension.
    Y: array
        The coordinates in the y dimension.
    Z: array
        The values of the function.
    elev: float
        The elevation angle in degrees rotates the camera above the plane
        pierced by the vertical axis, with a positive angle corresponding
        to a location above that plane.
    azim: float
        The azimuthal angle in degrees rotates the camera about the
        vertical axis, with a positive angle corresponding to a
        right-handed rotation.
    roll: float
        The roll angle in degrees rotates the camera about the viewing axis.

    Returns
    -------
        The matplotlib.axis object of the selected subplot.
    """
    my_cmap = plt.get_cmap('Spectral')
    ax = fig.add_subplot(nb_x, nb_y, n_idx+1, projection='3d')
    ax.view_init(elev, azim, roll)
    surf = ax.plot_trisurf(X, Y, Z, cmap = my_cmap, linewidth=0, antialiased=False)
    fig.colorbar(surf, ax = ax, shrink = 0.5, aspect = 5)
    return ax


def display(iter_nb, dt, title, subtitle, output_folder):
    """
    Plot the advected function for different time steps.

    Parameters
    ----------
    iter_nb: int
        Number of total time iteration of the simulation.
    dt: float
        The time step.
    title: string
        The name of the file where the figure is saved.
    subtitle: string
        The written name of the figure.
    """
    fig = plt.figure(figsize=(27,15))
    fig.suptitle(subtitle, fontsize=16)

    Idx_list = [0]+[min(i*(iter_nb//9) + 1,iter_nb) for i in range(9)] + [iter_nb]
    for n_idx in range(9 + 2):
        idx = Idx_list[n_idx]
        namefile = os.path.join(output_folder, f"after_{idx}.txt")
        list_F, _, _, CoordX, CoordY = treatment(namefile)

        X = np.array(CoordX)
        Y = np.array(CoordY)
        Z = np.array(list_F)

        ax = set_surface(fig, 3, 5, n_idx, X, Y, Z, 90, -90, 0)

        if idx == 0 :
            plt.title(f"Initial condition\n (maximum = {max(Z):.3}) ")
        elif idx <iter_nb-1:
            plt.title(f"Advected {idx} time"+"s"*(idx!=1)+f" (dt = {dt:.3})\n (maximum = {max(Z):.3}) ")
        else:
            plt.title(f"Final state ({iter_nb} times)\n (maximum = {max(Z):.3}) ")

        set_axis(ax, "x", "y", "z", "equal")


    namefile_final = os.path.join(output_folder, f"after_{iter_nb}.txt")
    namefile_exact = os.path.join(output_folder, f"after_{iter_nb}_exact.txt")
    list_F_final, _, _, CoordX, CoordY = treatment(namefile_final)
    list_F_exact, _, _, CoordX, CoordY = treatment(namefile_exact)


    X = np.array(CoordX)
    Y = np.array(CoordY)
    Zf = np.array(list_F_final)
    Ze = np.array(list_F_exact)
    Z = Ze - Zf

    ax = set_surface(fig, 3, 5, 11, X, Y, Ze, 90, -90, 0)
    plt.title("Expected final state")
    set_axis(ax, "x", "y", "z", "equal")

    ax = set_surface(fig, 3, 5, 12, X, Y, Z, 90, -90, 0)
    plt.title(f"Difference (=exact-advected)\n max error = {max(abs(Z)):.3} ")
    set_axis(ax, "x", "y", "z", "equal")

    ax = set_surface(fig, 3, 5, 13, X, Y, Z, 45, -45, 0)
    plt.title(f"Difference (=exact-advected)\n max error = {max(abs(Z)):.3} ")
    set_axis(ax, "x", "y", "z", "on")


    ax = set_surface(fig, 3, 5, 14, X, Y, -Z, 45, -45, 0)
    plt.title(f"Difference (=-(exact-advected))\n max error = {max(abs(Z)):.3} ")
    set_axis(ax, "x", "y", "-z", "on")


    #plt.savefig(title + ".png")
    plt.show()




# Get the inputs -----------------------------------------------
executable, yaml_parameters, _, verbose = set_input(0, 1, 60, 120, 0.01, 0.8,  True, False)

executable_name = os.path.basename(executable)
if executable_name == "advection_ALL" :
    print("Choose an executable of one simulation, not advection_ALL which launches severals simulations.")
assert(executable_name != "advection_ALL")


# Execute the test file given as input in the command ----------
out = execute(executable, yaml_parameters, verbose)


# Display the curves --------------------------------------------
mapping, method, domain, simulation, name = get_simulation_config(executable)

rmin = yaml_parameters['r_min']
rmax = yaml_parameters['r_max']
Nr = yaml_parameters['r_size']
Nth = yaml_parameters['p_size']
dt = yaml_parameters['time_step']
T = yaml_parameters['final_time']
iter_nb = int(T/dt)
details = f"\n $NrxNt$ = {Nr}x{Nth}; [$rmin$,$rmax$] = [{rmin},{rmax}]."
output_folder = f'{mapping.replace(" ","_")}_{domain}-{method.replace(" ","_")}-{simulation}_output'
display(iter_nb, dt, executable_name, name, output_folder)




