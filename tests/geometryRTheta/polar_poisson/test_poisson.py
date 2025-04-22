#!/usr/bin/python3
""" A file used for testing polarpoissonfemsolver.cpp. It provides 2 input files of different sizes and uses the results to calculate the order of convergence.
Cubic splines should produce 4th order convergence. As few points are used the convergence order will not yet have converged so we simply verify that order > 3.5.

Inputs: the executable associated with the file polarpoissonfemsolver.cpp.
"""
import subprocess
import sys
import pathlib

import numpy as np

executable = sys.argv[1]
exe_path=pathlib.PurePath(executable)
test_case=exe_path.stem.removeprefix('polar_poisson_convergence_')
input_file=test_case+".yaml"

with open(input_file, "w", encoding="utf-8") as f:
    print("SplineMesh:", file=f)
    print("  r_ncells: 16", file=f)
    print("  theta_ncells: 16", file=f)

with subprocess.Popen([executable, input_file], stdout=subprocess.PIPE) as p:
    out, err = p.communicate()
out = out.decode('ascii').strip()

assert p.returncode == 0
print(out)
if err:
    err = err.decode('ascii')
    print(err)

out_lines = out.split('\n')
error_32 = [float(l.split(' ')[3]) for l in out_lines if "Max error :" in l][0]

with open(input_file, "w", encoding="utf-8") as f:
    print("SplineMesh:", file=f)
    print("  r_ncells: 32", file=f)
    print("  theta_ncells: 32", file=f)

with subprocess.Popen([executable, input_file], stdout=subprocess.PIPE) as p:
    out, err = p.communicate()
out = out.decode('ascii').strip()

assert p.returncode == 0
print(out)
if err:
    err = err.decode('ascii')
    print(err)

out_lines = out.split('\n')
error_64 = [float(l.split(' ')[3]) for l in out_lines if "Max error :" in l][0]

order = np.log(error_32/error_64) / np.log(2)

print("Measured order : ", order)

assert 3.5 < order
