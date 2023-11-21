#!/usr/bin/python3
""" 
A file used for testing polarpoissonfemsolver.cpp.
It provides 2 input files of different sizes and uses the results 
to calculate the order of convergence.
Cubic splines should produce 4th order convergence. 
As few points are used the convergence order will not 
yet have converged so we simply verify that order > 3.5.

Inputs: the executable associated with the file polarpoissonfemsolver.cpp.
"""
import subprocess
import sys

import numpy as np

executable = sys.argv[1]

with open("poisson.yaml", "w", encoding="utf-8") as f:
    print("Mesh:", file=f)
    print("  r_size: 32", file=f)
    print("  p_size: 64", file=f)

with subprocess.Popen([executable, "poisson.yaml"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as p:
    out, err = p.communicate()
    assert p.returncode == 0

print(out)
if err:
    print(err)
    assert False

out_lines = out.split('\n')
error_64 = [float(l.split(' ')[4]) for l in out_lines if "Max error function : " in l][0]
error_dx_64 = [float(l.split(' ')[4]) for l in out_lines if "Max error Ex : " in l][0]
error_dy_64 = [float(l.split(' ')[4]) for l in out_lines if "Max error Ey : " in l][0]


with open("poisson.yaml", "w", encoding="utf-8") as f:
    print("Mesh:", file=f)
    print("  r_size: 64", file=f)
    print("  p_size: 128", file=f)

with subprocess.Popen([executable, "poisson.yaml"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as p:
    out, err = p.communicate()
    assert p.returncode == 0

print(out)
if err:
    print(err)
    assert False

out_lines = out.split('\n')
error_128 = [float(l.split(' ')[4]) for l in out_lines if "Max error function : " in l][0]
error_dx_128 = [float(l.split(' ')[4]) for l in out_lines if "Max error Ex : " in l][0]
error_dy_128 = [float(l.split(' ')[4]) for l in out_lines if "Max error Ey : " in l][0]


order = np.log(error_64/error_128) / np.log(2)
order_dx = np.log(error_dx_64/error_dx_128) / np.log(2)
order_dy = np.log(error_dy_64/error_dy_128) / np.log(2)

print("Measured order of the electric potential : ", order)
print("Measured order of the electric field in x : ", order_dx)
print("Measured order of the electric field in y : ", order_dy)


assert 3.5 < order
assert 2.5 < order_dx
assert 2.5 < order_dy
