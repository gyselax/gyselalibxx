#!/usr/bin/python3
import subprocess
import sys

import numpy as np

executable = sys.argv[1]

with open("poisson.yaml", "w") as f:
    print("Mesh:", file=f)
    print("  r_size: 64", file=f)
    print("  p_size: 64", file=f)

p = subprocess.Popen([executable, "poisson.yaml"], stdout=subprocess.PIPE)
out, err = p.communicate()
out = out.decode('ascii').strip()

assert p.returncode == 0
print(out)
if err:
    err = err.decode('ascii')
    print(err)
    assert False

out_lines = out.split('\n')
error_64 = float(out_lines[-1].split(' ')[-1])

with open("poisson.yaml", "w") as f:
    print("Mesh:", file=f)
    print("  r_size: 128", file=f)
    print("  p_size: 128", file=f)

p = subprocess.Popen([executable, "poisson.yaml"], stdout=subprocess.PIPE)
out, err = p.communicate()
out = out.decode('ascii').strip()

assert p.returncode == 0
print(out)
if err:
    err = err.decode('ascii')
    print(err)
    assert False

out_lines = out.split('\n')
error_128 = float(out_lines[-1].split(' ')[-1])

order = np.log(error_64/error_128) / np.log(2)

print("Measured order : ", order)

assert 3.5 < order
