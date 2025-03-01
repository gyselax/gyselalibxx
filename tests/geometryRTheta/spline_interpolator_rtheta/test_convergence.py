# Command to launch the test :
# python3 test_convergence.py ../../../build/tests/geometryRTheta/2d_spline_interpolator/2d_spline_interpolator_tests

"""
Compute the space convergence order between two grid sizes
of the 2D polar spline interpolator.

Parameters
----------
N0 : int
    Size of the first grid: N0xN0. 
    Possibility for the user to change it by adding the two 
    new sizes in the command line.

N1 : int
    Size of the second grid: N1xN1.
    Possibility for the user to change it by adding the two 
    new sizes in the command line.
   
Returns
-------
    List of convergence order for each test functions. 
"""


import argparse
import numpy as np

from polar_interpolation_functions import launch_executable


# Inputs --------------------------------------------------------
parser = argparse.ArgumentParser(description='Test the convergence order in space between two grid sizes.')
parser.add_argument('executable', metavar='executable', type=str, nargs= 1, help='path to be executable.')
parser.add_argument('N0', metavar='N0', type=int, nargs='?', default = 64, help='first grid size: N0xN0. By default, N0 = 64.')
parser.add_argument('N1', metavar='N1', type=int, nargs='?', default = 128, help='second grid size: N1xN1. By default, N1 = 128.')
args = parser.parse_args()

N0 = args.N0
N1 = args.N1

print("")


# Execute the test file given as input in the command ----------
executable = args.executable[0]

NN = [N0, N1]

[total_Errors0, total_Errors1], fct_names, degree = launch_executable(executable, NN)
nb_fct = len(fct_names)


# Compute the order ---------------------------------------------
nb_fct = len(fct_names)

print("Measured order : ")
Orders = np.zeros(nb_fct)
for i in range(nb_fct):
    Orders[i] = np.log(total_Errors0[i]/total_Errors1[i]) / np.log(N1/N0)
    print("  ", fct_names[i] +": ", Orders[i])

    # Order 4 expected:
    if degree +1 -0.25 > Orders[i]:
        if (total_Errors0[i] <1e-14 or total_Errors1[i] <1e-14):
            print("      WARNING: Errors very small: saturation with machine error.")
        else:
            print("      WARNING: The grid may not be refined enough.\n               Recommended to take as input N > 60 (for NxN grid).")
    assert (degree +1 -0.25 < Orders[i] or (total_Errors0[i] <1e-14 or total_Errors1[i] <1e-14))

print("Average measured order : ", sum(Orders)/len(Orders))


