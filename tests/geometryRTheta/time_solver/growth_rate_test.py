#!/usr/bin/python3

# SPDX-License-Identifier: MIT

"""
Inputs: the executable associated with the file name.cpp.
"""
# import subprocess
# import sys
# import pathlib

from argparse import ArgumentParser
from pathlib import Path

import numpy as np

from gysdata import DiskStore

# executable = sys.argv[1]
# exe_path = pathlib.PurePath(executable)

# test_case = exe_path.stem.removeprefix('polar_poisson_convergence_')
# path = ""
# input_file = path + "params.yaml"




# def check_L2_norms(folder):
if __name__ == '__main__':
    path_data_structure = Path('data_structure_RTheta.yaml')
    ds = DiskStore(folder, data_structure=path_data_structure)

    # Get initial data
    rho_eq = np.array(ds['density_eq'])
    phi_eq = np.array(ds['electrical_potential_eq'])

    jacobian = np.array(ds["jacobian"])

    # dt = float(ds["delta_t"])
    T = float(ds["final_T"])
    # iter_nb = int( T / dt ) +1
    # tstep_diag = int(ds["time_step_diag"])

    # Nr = int(ds["r_size"])
    # Nt = int(ds["theta_size"])
    omega_Im = float(ds["slope"])

    # Get the data at each time step
    Time = np.array(ds["density"].coords["time"])
    rho = np.array(ds['density'])
    phi = np.array(ds['electrical_potential'])

    # Compute norms
    L2norms_rho = np.linalg.norm((rho - rho_eq[None,:,:])*abs(jacobian[None,:,:]), 2, axis=(1,2))
    L2norms_phi = np.linalg.norm((phi - phi_eq[None,:,:])*abs(jacobian[None,:,:]), 2, axis=(1,2))

    slope_rho = np.polyfit(Time, L2norms_rho, 1)[0]
    slope_phi = np.polyfit(Time, L2norms_phi, 1)[0]

    assert abs(slope_rho - omega_Im) < 1e-3
    assert abs(slope_phi - omega_Im) < 1e-3