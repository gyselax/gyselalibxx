# Command to launch the test :
# python3 display_errors.py ../../../build/tests/geometryRTheta/2d_spline_interpolator/2d_spline_interpolator_tests

"""
Display the maximum absolute errors of each test functions
made by the 2D polar spline interpolator for several grid sizes
and compute the slope of each curves. 

Parameters
----------
executable: string
    Path to the executable.

Returns
-------
    : matplotlib objet
    The curves are printed. 
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

from polar_interpolation_functions import launch_executable

# Inputs --------------------------------------------------------
NN = range(16, 128+1)


# Execute the test files given as input in the command ----------
executable = sys.argv[1]

total_Errors, fct_names, degree = launch_executable(executable, NN)
nb_fct = len(fct_names)


# Compute the order ---------------------------------------------
### Sort by function test:
Errors_by_funct = [[total_Errors[j][i] for j in range(len(NN))] for i in range(nb_fct)]


# Plot the errors curves ----------------------------------------
plt.figure(figsize = (8,6))
for i in range(nb_fct):
    max_index = 0
    while Errors_by_funct[i][max_index] > 1e-14 and max_index < len(NN)-1:
        max_index += 1
    slope, intercept = np.polyfit(-np.log(NN[:max_index]),np.log(Errors_by_funct[i][:max_index]),1)
    plt.loglog(NN, Errors_by_funct[i], "+--", label=f"{fct_names[i]} (slope = {round(slope,2)})", linewidth = 0.5)

    # Order d+1 expected:
    assert (slope > degree + 1 - 0.25)

plt.legend()
plt.xlabel("N for NxN grid (log)")
plt.ylabel("Errors (log)")
plt.grid()
plt.title(f"Space errors of interpolator (r, theta) on bsplines of degree {degree}")
plt.show()



