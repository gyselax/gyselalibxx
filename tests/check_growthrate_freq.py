#!/bin/env python3

# SPDX-License-Identifier: MIT

""" File which tests whether the specified results of a landau damping test case
follow the gradient predicted by the theory
"""

import sys

from argparse import ArgumentParser
from pathlib import Path

import diag_efield as Ediag
from gysdata import DiskStore

if __name__ == '__main__':
    parser = ArgumentParser(description="Plot frequency.")
    parser.add_argument('data_dir',
                        action='store',
                        nargs='?',
                        default=Path.cwd(),
                        type=Path,
                        help='location of the results')
    parser.add_argument('-g', '--growthrate',
                        action='store',
                        default=None,
                        type=float,
                        help='theoretical growth rate')
    parser.add_argument('-f', '--frequency',
                        action='store',
                        default=None,
                        type=float,
                        help='theoretical frequency')
    args = parser.parse_args()

    # Load data
    ds = DiskStore(args.data_dir)
    epot = ds['electrostatic_potential']
    epot = epot.isel(x=int(epot.sizes["x"] / 2))

    # Compute and plot the growth (or damping) rate
    (fitted_values, fit,
     growthrate_computed, valid_growthrate) = Ediag.compute_growthrate(epot,
                                                                       args.growthrate)
    growthrate_outfile = Path(f'growthrate_t{epot.coords["time"].values[0]}'
                              f'to{epot.coords["time"].values[-1]}.png')
    Ediag.plot_growthrate(growthrate_outfile, abs(epot), fitted_values, fit,
                          growthrate_computed, args.growthrate)

    # Compute and plot the growth (or damping) rate
    (frequencies, max_frequency, valid_freq) = Ediag.compute_frequency(epot,
                                                                       args.frequency)
    frequency_outfile = Path(f'frequency_t{epot.coords["time"].values[0]}'
                             f'to{epot.coords["time"].values[-1]}.png')
    Ediag.plot_frequency(frequency_outfile, frequencies, max_frequency, args.frequency)

    if not valid_growthrate and not valid_freq:
        sys.exit(1)
