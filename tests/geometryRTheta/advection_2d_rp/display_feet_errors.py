# Command to launch the test :
# python3 display_feet_errors.py ../../../build/tests/geometryRTheta/advection_2d_rp/advection_2d_rp_tests

"""
Display the characteristic feet and their errors computed at the last time step of the advection
simulation.

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

from advection_functions import set_input, execute, params_file, treatment_feet, distance


# Definition of functions ---------------------------------------
# --- Compute errors
def compute_max_distance_error(Nr, Nt, computed_Fx, computed_Fy,exact_Fx, exact_Fy):
    """
    Compute the maximum distance betweent the exact feet and the computed feet.

    Parameters
    ----------
    Nr: int
        The number of cells in the r dimension.
    Nt: int
        The number of cells in the theta dimension.
    computed_Fx: list
        The computed characteristic feet in x.
    computed_Fy: list
        The computed characteristic feet in y.
    exact_Fx: list
        The exact characteristic feet in x.
    exact_Fy: list
        The exact characteristic feet in y.

    Returns
    -------
        The maximum distance between the computed feet and the exact feet.
    """
    # Troncate the list to avoid the errors on the edge due to the evaluator used:
    size = Nr*Nt // 2
    D = distance(exact_Fx[:size], exact_Fy[:size], computed_Fx[:size], computed_Fy[:size])
    err_d = np.sqrt(sum(D*D))
    return err_d

def compute_distance_errors(Nr, Nt, computed_Fx, computed_Fy, exact_Fx, exact_Fy):
    """
    Compute the distances between the exact feet and the computed feet.

    Parameters
    ----------
    Nr: int
        The number of cells in the r dimension.
    Nt: int
        The number of cells in the theta dimension.
    computed_Fx: list
        The computed characteristic feet in x.
    computed_Fy: list
        The computed characteristic feet in y.
    exact_Fx: list
        The exact characteristic feet in x.
    exact_Fy: list
        The exact characteristic feet in y.

    Returns
    -------
        A list of the distances between the computed feet and the exact feet.
    """
    distances = distance(exact_Fx, exact_Fy, computed_Fx, computed_Fy)
    err_d = np.reshape(distances, (Nr,Nt))
    return err_d


# --- Set plots
def set_plot_slope(idx, Values_list, title):
    """
    Plot the maximum distance errors and the slope of the errors.

    Parameters
    ----------
    idx: int
        The index of the subplot.
    Values_list: list
        A list of the errors.
    title: string
        The title of the subplot.
    """
    plt.subplot(idx)

    # Compute the slope while the values are not null:
    if Values_list != [0.]*len(Values_list):
        index = 0
        while index < len(Values_list) and Values_list[index] != 0 :
            index+=1
        slope = np.polyfit(np.log(DT[:index]),np.log(Values_list[:index]),1)
        plt.loglog(DT, Values_list, "+--", linewidth = 0.5, label=f"Max errors (slope = {round(slope[0],2)})" )

    # Do not compute the slope if the values are null:
    else :
        plt.loglog(DT, Values_list, "+--", linewidth = 0.5, label="Max errors (null errors)" )

    plt.legend()
    plt.xlabel("dt (log)")
    plt.ylabel("Errors (log)")
    plt.grid()
    plt.xticks(DT)
    plt.title(title)


def set_plot_imshow(index, Values_list, title):
    """
    Plot the distance errors on the domain.

    Parameters
    ----------
    idx: int
        The index of the subplot.
    Values_list: list
        A list of the errors.
    title: string
        The title of the subplot.
    """
    plt.subplot(index)
    plt.imshow(Values_list)
    plt.colorbar()
    plt.xlabel("theta")
    plt.ylabel("r")
    plt.title(title)



# Get the inputs ------------------------------------------------
executable, rmin, rmax, Nr, Nt, dt0, T, curves, feet, _ = set_input(0, 1, 40, 80, 0.1, 0.8,  False, True)

executable_name = os.path.basename(executable)
if executable_name == "advection_ALL" :
    print("Choose an executable of one simulation, not advection_ALL which launches severals simulations.")
assert(executable_name != "advection_ALL")

max_pow = 9
DT = [dt0 * 2**(-i) for i in range(0,max_pow)]



ERR_d = []
ERR_d_init_end = []
L_exact_Fr, L_exact_Fp, L_computed_Fr, L_computed_Fp  = [], [], [], []
for dt in DT:
    # Execute the test file given as input in the command -------
    execute(executable, rmin, rmax, Nr, Nt, dt, T, curves, True)

    namefile1 = "output/feet_computed_-1.txt"
    namefile2 = "output/feet_exact_-1.txt"

    _, _, Cr, Cp, _, _, computed_Fr, computed_Fp, computed_Fx, computed_Fy  = treatment_feet(namefile1)
    _, _, _, _, _, _,  exact_Fr, exact_Fp, exact_Fx, exact_Fy = treatment_feet(namefile2)

    # --- Plot type 1 : the convergence order
    ERR_d += [compute_max_distance_error(Nr+3, Nt, computed_Fx, computed_Fy, exact_Fx, exact_Fy)]

    # --- Plot type 2 : errors on the domain
    if dt in (DT[0], DT[-1]):
        ERR_d_init_end += [compute_distance_errors(Nr+3, Nt, computed_Fx, computed_Fy, exact_Fx, exact_Fy)]

    # --- Plot type 3 : the computed and the exact feet
    L_exact_Fr +=[exact_Fr]
    L_exact_Fp +=[exact_Fp]
    L_computed_Fr +=[computed_Fr]
    L_computed_Fp +=[computed_Fp]




# Put "False" the savings of files for the next launch ----------
params_file(rmin, rmax, Nr, Nt, dt0, T)


# Name of the figure --------------------------------------------
executable_name = os.path.basename(executable)
name = executable_name.lower().split('__')
fct_names = name[2].lower().replace("_", " ")
fct_names += " with " + name[1].lower().replace("_", " ")
fct_names += " on " + ' '.join(name[0].split("_")[1:]).lower()
fct_names = fct_names[0].upper() + fct_names[1:]

# Plot the feet -------------------------------------------------
fig = plt.figure(figsize = (20,20))
fig.suptitle(fct_names, fontsize=16)


set_plot_slope(341, ERR_d,  f"Time errors on distance between \n exact and computed ({Nr}x{Nt})")
set_plot_imshow(345, ERR_d_init_end[0],  f"Distance errors for dt = {DT[0]}")
set_plot_imshow(349, ERR_d_init_end[-1],  f"Distance errors for dt = {DT[-1]}")


for i in range(max_pow):
    ax =plt.subplot(3, 4, i+1 +i//3 +1)
    plt.plot(Cr, Cp, "+", color="silver", label="grid")
    plt.plot(L_exact_Fr[i], L_exact_Fp[i], "+", label="exact")
    plt.plot(L_computed_Fr[i], L_computed_Fp[i], "+", label="computed")

    plt.xlabel("r")
    plt.ylabel("theta")
    plt.legend()
    plt.title(f"Exact and computed feet \n for dt = {DT[i]} and NrxNt = {Nr}x{Nt}", fontdict={'fontsize': 10 })
    ax.set_adjustable("box")

plt.show()





