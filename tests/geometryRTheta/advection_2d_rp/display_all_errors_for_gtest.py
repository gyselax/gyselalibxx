# Command to launch the test :
# python3 display_all_error_for_gtest.py ../../../build/tests/geometryRTheta/advection_2d_rp/advection_ALL

"""
Compute the space convergence order between two grid sizes
of the 2D polar spline interpolator.

Parameters
----------
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


Returns
-------
    The convergence orders for each test functions.
"""
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt

from advection_functions import set_input, execute, extract_simulation_results



# Get the inputs -----------------------------------------------
executable, yaml_parameters, plot_results, verbose = set_input(0, 1, 32, 64, 0.1, 0.4)

if plot_results:
    list_pow = [0, 1, 2, 3]
else :
    list_pow = [0, 2]

yaml_configs = [deepcopy(yaml_parameters) for _ in list_pow]

Nr0 = yaml_parameters['SplineMesh']['r_ncells']
Nt0 = yaml_parameters['SplineMesh']['p_ncells']
dt0 = yaml_parameters['Time']['time_step']
for i, p in enumerate(list_pow):
    yaml_configs[i]['Time']['time_step'] = dt0 / 2**p
    yaml_configs[i]['SplineMesh'].update({'r_ncells': Nr0 * 2**p,
                                          'p_ncells': Nt0 * 2**p})
Order = [2**-i for i in list_pow]


# Expected convergence orders --------------------------
order_expected = {'euler': {
                        'circular physical': {'Translation':3, 'Rotation':1, 'Decentred rotation':1},
                        'czarny physical': {'Translation':3, 'Rotation':1, 'Decentred rotation':1},
                        'czarny pseudo cartesian': {'Translation':1, 'Rotation':1, 'Decentred rotation':1},
                        'discrete pseudo cartesian': {'Translation':1, 'Rotation':1, 'Decentred rotation':1}
                      },
                  'crank nicolson': {
                        'circular physical': {'Translation':3, 'Rotation':2, 'Decentred rotation':2},
                        'czarny physical': {'Translation':3, 'Rotation':2, 'Decentred rotation':2},
                        'czarny pseudo cartesian': {'Translation':3, 'Rotation':2, 'Decentred rotation':2},
                        'discrete pseudo cartesian': {'Translation':3, 'Rotation':2, 'Decentred rotation':2}
                    },
                  'rk3': {
                        'circular physical': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'czarny physical': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'czarny pseudo cartesian': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'discrete pseudo cartesian': {'Translation':3, 'Rotation':3, 'Decentred rotation':3}
                    },
                  'rk4': {
                        'circular physical': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'czarny physical': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'czarny pseudo cartesian': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'discrete pseudo cartesian': {'Translation':3, 'Rotation':3, 'Decentred rotation':3}
                        }
                  }

# Run the simulations -------------------------------------------
execution_output = [execute(executable, params, verbose) for params in yaml_configs]

# Get the results of the simulation
simulation_results = extract_simulation_results(execution_output)
nb_fct = len(simulation_results)

# Check the convergence order -----------------------------------
for (method, mapping, problem_type), result in simulation_results.items():
    Coeffs = np.polyfit(np.log(Order), np.log(result.l_inf_error),1)
    slope = Coeffs[0]

    theoretical_order = order_expected[method][mapping][problem_type]
    print(f">{result.name} :")
    print(f"      order = {round(slope,3)} ({theoretical_order} expected).")

    # Order expected:
    if (slope < theoretical_order - 0.25):
        print(f">{result.name} :")
        print(f" Wrong convergence order: got {slope}, expected {theoretical_order - 0.25}.")
    assert (slope > theoretical_order - 0.25)

# Plot the errors -----------------------------------------------
if plot_results:
    fig = plt.figure(figsize = (16,9))
    fig.suptitle(f"Errors of advection operator (starting with dt = {dt0}$^*$ and NrxNt = {Nr0}x{Nt0})\n ($^*$dt = {dt0/10} for Euler method)", fontsize=16)

    def set_subplot(index, method, mapping):
        """
        Plot the errors and their slope for the three
        test simulations: translated Gaussianne;
        rotated Gaussianne and Edoardo's test case.

        Parameters
        ----------
        index: int
            Index of the subplot.
        method : str
            The method used.
        mapping : str
            The mapping used.
        """
        axe = fig.add_subplot(4, 4, index)
        plt.tight_layout(pad=0.4, w_pad=0.2, h_pad=-0.2)


        for simu in ('Translation', 'Rotation', 'Decentred rotation'):
            try:
                _, linf_errors, _ = simulation_results[(method, mapping, simu)]
            except KeyError:
                print(f"Simulation {method}, {mapping}, {simu} was not run")
                continue
            slope = np.polyfit(np.log(Order),np.log(linf_errors),1)
            axe.loglog(Order, linf_errors, "+--", label=f"{simu} (slope = {round(slope[0],2)})", linewidth = 0.7)

        axe.legend()
        axe.grid()

        axe.set_xticks(Order)
        axe.set_ylabel("Errors (log)")
        axe.set_adjustable("box")

        axe.set_title(f'{method.capitalise()} on {mapping}', fontsize=12)

    index = 1
    for method, method_types in order_expected.items():
        for mapping in method_types:
            set_subplot(index, method, mapping)
            index += 1


    plt.show()



