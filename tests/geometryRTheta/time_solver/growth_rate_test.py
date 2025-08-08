#!/usr/bin/python3

# SPDX-License-Identifier: MIT

"""
Inputs: the executable associated with the file name.cpp.
"""

from argparse import ArgumentParser
from pathlib import Path
import h5py

import numpy as np

from gysdata import DiskStore


if __name__ == '__main__':
    parser = ArgumentParser(description="Check growth rate.")
    parser.add_argument('data_dir',
                        action='store',
                        nargs='?',
                        default=Path.cwd(),
                        type=Path,
                        help='Location of the results.')
    args = parser.parse_args()

    path_data_structure = Path('data_structure_RTheta.yaml')
    # folder = args.data_dir.joinpath("/output/")
    folder = Path.cwd()/"output"
    print("path_data_structure =", path_data_structure)
    print("Path.cwd() =", Path.cwd())
    print("args.data_dir =", args.data_dir)
    print("folder =", folder)
    data = h5py.File(folder/"GYSELALIBXX_initstate.h5")
    print ("keys:", list(data.keys()))

    ds = DiskStore(folder, data_structure=path_data_structure)
    print(ds)
    print(ds.keys())

    # Get initial data
    rho_eq = np.array(ds['density_eq'])
    phi_eq = np.array(ds['electrical_potential_eq'])

    jacobian = np.array(ds["jacobian"])

    T = float(ds["final_T"])

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

    print("abs(slope_rho - omega_Im)", abs(slope_rho - omega_Im))
    print("abs(slope_phi - omega_Im)", abs(slope_phi - omega_Im))

    assert abs(slope_rho - omega_Im) < 1e-3
    assert abs(slope_phi - omega_Im) < 1e-3
