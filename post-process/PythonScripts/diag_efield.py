# SPDX-License-Identifier: MIT

'''
Diagnostics associated to the electrostatic potential and electric field
'''

# pylint: disable=import-error
# pylint: disable=no-member

import sys

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import math_utils as mathut


def compute_growthrate(epot: xr.DataArray,
                       theoretical_growthrate=None):
    '''
    Compute the growth rate (or damping rate) from the time
    evolution of abs(Phi)
    '''

    valid = True
    # select the log of the value of local extrema (maxima in abs value)
    abs_epot = abs(epot)
    epot_extrema = abs_epot.where(
        np.logical_and(
            (abs_epot - abs_epot.shift(time=-1)) > 0,
            (abs_epot - abs_epot.shift(time=1)) > 0).compute(),
        drop=True)

    # find a linear fit for the log of the extrema
    # skip points at the end to ensure that the correlation coef
    # is close to 1 (10^-3 error max)
    fitted_values = np.log(epot_extrema[3:])
    while (np.square(np.corrcoef(fitted_values.coords['time'],
                                 fitted_values)) < 0.999).any():
        fitted_values = fitted_values[:-1]
    fit = fitted_values.polyfit("time", 1)
    fitted_values = np.exp(fitted_values)
    growthrate = fit['polyfit_coefficients'].values[0]
    fit = np.exp(fit['polyfit_coefficients'].values[0] *
                 abs_epot.coords['time']
                 + fit['polyfit_coefficients'].values[1])

    print(f'Computed growthrate: {growthrate:0.4f}')
    if theoretical_growthrate is not None:
        growthrate_relerror = abs(growthrate - theoretical_growthrate)
        growthrate_relerror = growthrate_relerror/abs(theoretical_growthrate)
        print(f'Expected ~{theoretical_growthrate:0.3f}:'
              f' error: {growthrate_relerror*100:0.1f}%', file=sys.stderr)
        if growthrate_relerror > 0.05:
            valid = False

    return (fitted_values, fit, growthrate, valid)


def plot_growthrate(filename, epot, fitted_values, fit, growthrate,
                    theoretical_growthrate=None):
    '''
    Plot the time evolution of abs(Phi)
    '''

    plt.figure()
    plt.semilogy(epot.coords["time"], epot)
    plt.semilogy(fit.coords["time"], fit, 'r', lw=2)
    plt.semilogy(fitted_values.coords["time"],
                 fitted_values, 'o', markersize=6)
    plt.ylabel('abs(Phi(middle_box))', fontsize=16)
    plt.xlabel('time', fontsize=16)
    if theoretical_growthrate is not None:
        plt.title(f'growthrate: {growthrate:0.3f} (theory: {theoretical_growthrate:0.3f})',
                  fontsize=16)
    else:
        plt.title(f'growthrate: {growthrate:0.3f}', fontsize=16)
    print("Save figure in " / filename)
    plt.savefig(filename)


def compute_frequency(epot: xr.DataArray,
                      theoretical_frequency=None):
    '''
    Compute the frequency from the time evolution of abs(Phi)
    '''

    valid = True
    frequencies = abs(mathut.fft1d(epot))
    max_frequency = abs(frequencies.idxmax('time_freq').values[()])
    print(f'Computed frequency: {max_frequency:0.4f}')

    if theoretical_frequency is not None:
        frequency_relerror = abs(max_frequency - theoretical_frequency)
        frequency_relerror = frequency_relerror / theoretical_frequency
        print(f'Expected ~{theoretical_frequency:0.3f}:'
              f' error: {frequency_relerror*100:0.1f}%', file=sys.stderr)
        if frequency_relerror > 0.05:
            valid = False

    return (frequencies, max_frequency, valid)


def plot_frequency(filename, frequencies, max_frequency,
                   theoretical_frequency=None):
    '''
    Plot the frequency from the time evolution of abs(Phi)
    '''

    plt.figure()
    y_coords = [min(abs(frequencies)), max(abs(frequencies))]
    plt.plot(frequencies.coords['time_freq'], frequencies)
    plt.plot([max_frequency, max_frequency], y_coords, 'r-', lw=3,
             label=f'Simulation: {max_frequency:0.4f}')
    if theoretical_frequency is not None:
        plt.plot([theoretical_frequency, theoretical_frequency],
                 y_coords, 'g--', lw=3,
                 label=f'Theory: {theoretical_frequency:0.4f}')
    plt.title(f'omega: {max_frequency:0.4f}', fontsize=16)
    plt.ylabel('FFT(Phi(t))', fontsize=16)
    plt.xlabel('Freq', fontsize=16)
    plt.legend()

    print("Save figure in " / filename)
    plt.savefig(filename)


def plot_epot_atonetime(epot: xr.DataArray,
                        itime=-1,
                        output_filename=None):
    '''
    Plot the electrostatic potential at one specific time
    '''

    nb_time = epot.sizes["time"]
    assert(itime < nb_time)
    itime_array = [0, itime, -1]
    time_diag = epot.coords["time"].values[itime]

    line_type = ['-b', '-r', '-g']

    plt.figure()
    _, axs = plt.subplots(1, 3, figsize=(18, 6))

    # plot 1: E(t,x)
    ix1 = 0
    axs[ix1].pcolormesh(epot.coords['x'],
                        epot.coords['time'],
                        epot,
                        shading='auto',
                        cmap='jet')
    axs[ix1].set_title('electrostatic potential')
    axs[ix1].set_xlabel(f"${epot.coords['x'].name}$", fontsize=12)
    axs[ix1].set_ylabel(f"${epot.coords['time'].name}$", fontsize=12)

    # plot 2: time evolution of E(x=Nx/2)
    ix2 = 1
    nb_x = epot.sizes["x"]
    axs[ix2].semilogy(epot.coords["time"], np.abs(epot.isel(x=nb_x//2)))
    axs[ix2].set_xlabel(f"${epot.coords['time'].name}$", fontsize=12)
    axs[ix2].set_ylabel('abs(Phi(middle_box))', fontsize=12)

    # plot 3: E(x) for 3 times
    ix3 = 2
    for itime_index, itime_val in enumerate(itime_array):
        axs[ix3].plot(epot.coords['x'],
                      epot.isel(time=itime_val),
                      line_type[itime_index],
                      label=f't= {epot.coords["time"].values[itime]}')
    axs[ix3].set_xlabel(f"${epot.coords['x'].name}$", fontsize=12)
    axs[ix3].set_ylabel('epot', fontsize=12)
    axs[ix3].legend()

    if output_filename is None:
        output_filename = Path(f'epot_t{time_diag}.png')
    print("Save figure in " / output_filename)
    plt.savefig(output_filename)

