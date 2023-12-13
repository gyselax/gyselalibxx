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


import numpy as np
import matplotlib.pyplot as plt

from advection_functions import set_input, execute, take_errors_from_out, get_full_fct_names, get_fct_names, get_fct_name_key



# Get the inputs -----------------------------------------------
executable, rmin, rmax, Nr0, Nt0, dt0, T, curves, feet, if_plot = set_input(0, 1, 32, 64, 0.1, 0.4)

if if_plot:
    list_pow = [0, 1, 2, 3]
else :
    list_pow = [0, 2]

NNr = [Nr0 * 2**i for i in list_pow]
NNt = [Nt0 * 2**i for i in list_pow]
DT =  [dt0 / 2**i for i in list_pow]
Order = [2**(-i) for i in list_pow]



# Get the names of the test cases -------------------------------
out = execute(executable, rmin, rmax, 4, 4, T, T, False, False, False)
fct_names = get_full_fct_names(out)
short_fct_names = get_fct_names(out)
nb_fct = len(fct_names)

# Expected convergence orders --------------------------
order_expected = {'euler': {
                        'circular': {'Translation':3, 'Rotation':1, 'Decentred rotation':1},
                        'czarny_physical': {'Translation':3, 'Rotation':1, 'Decentred rotation':1},
                        'czarny_pseudo_cartesian': {'Translation':1, 'Rotation':1, 'Decentred rotation':1},
                        'discrete': {'Translation':1, 'Rotation':1, 'Decentred rotation':1}
                      },
                  'crank nicolson': {
                        'circular': {'Translation':3, 'Rotation':2, 'Decentred rotation':2},
                        'czarny_physical': {'Translation':3, 'Rotation':2, 'Decentred rotation':2},
                        'czarny_pseudo_cartesian': {'Translation':3, 'Rotation':2, 'Decentred rotation':2},
                        'discrete': {'Translation':3, 'Rotation':2, 'Decentred rotation':2}
                    },
                  'rk3': {
                        'circular': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'czarny_physical': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'czarny_pseudo_cartesian': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'discrete': {'Translation':3, 'Rotation':3, 'Decentred rotation':3}
                    },
                  'rk4': {
                        'circular': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'czarny_physical': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'czarny_pseudo_cartesian': {'Translation':3, 'Rotation':3, 'Decentred rotation':3},
                        'discrete': {'Translation':3, 'Rotation':3, 'Decentred rotation':3}
                        }
                  }



# Get the errors ------------------------------------------------
total_Errors = [take_errors_from_out(execute(executable, rmin, rmax, NNr[i], NNt[i], DT[i], T, False, False, False)) for i in range(len(Order))]
Errors_by_fct = [[err[j] for err in total_Errors] for j in range(nb_fct)]


# Check the convergence order -----------------------------------
for index in range(nb_fct):
    Coeffs = np.polyfit(np.log(Order),np.log(Errors_by_fct[index]),1)
    slope = Coeffs[0]

    problem_type, time_interation_method, mapping = get_fct_name_key(fct_names[index])

    theoretical_order = order_expected[time_interation_method][mapping][problem_type]
    print(f">{fct_names[index]} : \n      order = {round(slope,3)} ({theoretical_order} expected).")

    # Order expected:
    if (slope < theoretical_order - 0.25):
        print(fct_names[index] + f" has not the right convergence order: got {slope}, expected {theoretical_order - 0.25}.")
    assert (slope > theoretical_order - 0.25)



# Plot the errors -----------------------------------------------
if if_plot:
    fig = plt.figure(figsize = (16,9))
    fig.suptitle(f"Errors of advection operator (starting with dt = {dt0}$^*$ and NrxNt = {Nr0}x{Nt0})\n ($^*$dt = {dt0/10} for Euler method)", fontsize=16)

    def set_subplot(index, i0, title):
        """
        Plot the errors and their slope for the three
        test simulations: translated Gaussianne;
        rotated Gaussianne and Edoardo's test case.

        Parameters
        ----------
        index: int
            Index of the subplot.
        i0: int
            Index corresponding to the translation test case.
        title : string
            The title of the method and the mapping used.
        """
        ax = fig.add_subplot(4, 4, index)
        plt.tight_layout(pad=0.4, w_pad=0.2, h_pad=-0.2)

        # Translation
        slope = np.polyfit(np.log(Order),np.log(Errors_by_fct[i0]),1)
        ax.loglog(Order, Errors_by_fct[i0], "+--", label=f"{short_fct_names[i0]} (slope = {round(slope[0],2)})", linewidth = 0.7)

        # Rotation
        slope = np.polyfit(np.log(Order),np.log(Errors_by_fct[i0+1]),1)
        ax.loglog(Order, Errors_by_fct[i0+1], "+--", label=f"{short_fct_names[i0+1]} (slope = {round(slope[0],2)})", linewidth = 0.7)

        # Decentred rotation
        slope = np.polyfit(np.log(Order),np.log(Errors_by_fct[i0+2]),1)
        ax.loglog(Order, Errors_by_fct[i0+2], "+--", label=f"{short_fct_names[i0+2]} (slope = {round(slope[0],2)})", linewidth = 0.7)

        ax.legend()
        ax.grid()

        ax.set_xticks(Order)
        ax.set_ylabel("Errors (log)")
        ax.set_adjustable("box")

        ax.set_title(title, fontsize=12)


    names_method = ["Euler", "Crank Nicolson", "RK3", "RK4"]
    names_mapping = ["Circular mapping and physical domain",
                     "Czarny mapping and physical domain",
                     "Czarny mapping and pseudo Cartesian domain",
                     "Discrete mapping and pseudo Cartesian domain"]

    for i in range(nb_fct//3):
        index = i + 1
        set_subplot(index, 3*i, names_mapping[i//4] + " \n " + names_method[i%4])


    plt.show()



