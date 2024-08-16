#!/usr/bin/python3
""" A file used for testing polarpoissonfemsolver.cpp. It provides 2 input files of different sizes and uses the results to calculate the order of convergence.
Cubic splines should produce 4th order convergence. As few points are used the convergence order will not yet have converged so we simply verify that order > 3.5.

Inputs: the executable associated with the file polarpoissonfemsolver.cpp.
"""
import subprocess
import sys

import numpy as np

executable = sys.argv[1]

with open("poisson.yaml", "w", encoding="utf-8") as f:
    print("SplineMesh:", file=f)
    print("  r_ncells: 64", file=f)
    print("  p_ncells: 64", file=f)

with subprocess.Popen([executable, "poisson.yaml"], stdout=subprocess.PIPE) as p:
    out, err = p.communicate()
out = out.decode('ascii').strip()

assert p.returncode == 0
print(out)
if err:
    err = err.decode('ascii')
    print(err)

out_lines = out.split('\n')
error_64 = [float(l.split(' ')[3]) for l in out_lines if "Max error :" in l][0]

with open("poisson.yaml", "w", encoding="utf-8") as f:
    print("SplineMesh:", file=f)
    print("  r_ncells: 128", file=f)
    print("  p_ncells: 128", file=f)

with subprocess.Popen([executable, "poisson.yaml"], stdout=subprocess.PIPE) as p:
    out, err = p.communicate()
out = out.decode('ascii').strip()

assert p.returncode == 0
print(out)
if err:
    err = err.decode('ascii')
    print(err)

out_lines = out.split('\n')
error_128 = [float(l.split(' ')[3]) for l in out_lines if "Max error :" in l][0]

order = np.log(error_64/error_128) / np.log(2)

print("Measured order : ", order)

assert 3.5 < order
