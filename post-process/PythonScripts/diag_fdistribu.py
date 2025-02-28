# SPDX-License-Identifier: MIT

'''
Diagnostics associated to the distribution function
'''

# pylint: disable=import-error

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


def compute_deltaf(fdistribu_eq: xr.DataArray, fdistribu: xr.DataArray):
    '''
    Compute delta f function = fdistribu - fdistribu_eq
    '''
    deltaf = (fdistribu - fdistribu_eq)
    deltaf.name = f'{fdistribu.name}-{fdistribu_eq.name}'

    return deltaf


def plot_fdistribu_atfourtimes(fdistribu: xr.DataArray, figurename):
    '''
    Plot the distribution function at four time along the simulation
    '''

    plt.figure()
    nb_time = fdistribu.sizes["time"]
    nb_x = fdistribu.sizes["x"]
    ix_array = np.linspace(nb_x//8, nb_x, 3, dtype=int, endpoint=False)
    itime_array = np.linspace(0, nb_time-1, 4, dtype=int, endpoint=True)
    style_true = ['b-', 'c-', 'r-', 'k-']
    _, axs = plt.subplots(1, len(ix_array), figsize=(30, 8))
    for iix, ix in enumerate(ix_array):
        str_x_iix = str(round(float(fdistribu.coords['x'][ix]), 3))
        for itime_index, itime_val in enumerate(itime_array):
            str_time = str(round(float(fdistribu.coords['time'][itime_val]), 3))
            axs[iix].plot(fdistribu.coords['v_x'],
                          fdistribu.isel(time=itime_val, x=ix),
                          style_true[itime_index],
                          label='t = '+str_time)
            subtitle = 'x = '+str_x_iix
            axs[iix].set_title(subtitle, fontsize=16)
            axs[iix].set_xlabel('v_x', fontsize=16)
            if iix == 0:
                axs[iix].set_ylabel('Dist. function', fontsize=16)
            axs[iix].legend()

    print("Save figure in " + figurename)
    plt.savefig(figurename)


def plot_fdistribu_atonetime(fdistribu: xr.DataArray,
                             itime=-1,
                             output_filename=None):
    '''
    Plot the distribution function at one specific time
    '''

    nb_time = fdistribu.sizes["time"]
    assert(itime < nb_time)
    itime_array = [0, itime, -1]
    time_diag = fdistribu.coords["time"].values[itime]

    nb_x = fdistribu.sizes["x"]
    ix = nb_x//2
    nb_species = fdistribu.sizes["species"]
    nb_vx = fdistribu.sizes["v_x"]
    ivx_array = np.linspace(nb_vx//8, nb_vx, 3, dtype=int, endpoint=False)
    line_type = ['-b', '-r', '-g']

    plt.figure()
    _, axs = plt.subplots(nb_species, 3, figsize=(18, 6*nb_species))
    for ispec, _ in enumerate(fdistribu.coords['species']):
        f_spec = fdistribu.isel(time=itime, species=ispec)
        f_spec_name = f'fdistribu {fdistribu.coords["species"].item()}'

        # plot 1: f(x,v_x) at t=time_diag
        ix1 = 3*ispec
        axs[ix1].pcolormesh(f_spec.coords['v_x'],
                            f_spec.coords['x'],
                            f_spec.values,
                            shading='auto',
                            cmap='jet')
        axs[ix1].set_title(f'{f_spec_name} at t = {time_diag}')
        axs[ix1].set_xlabel(f"${f_spec.coords['v_x'].name}$", fontsize=12)
        axs[ix1].set_ylabel(f"${f_spec.coords['x'].name}$", fontsize=12)

        # plot 2: f(v) for 3 times
        ix2 = 3*ispec+1
        for itime_index, itime_val in enumerate(itime_array):
            axs[ix2].plot(fdistribu.coords['v_x'],
                          fdistribu.isel(time=itime_val, species=ispec, x=ix),
                          line_type[itime_index],
                          label=f't= {fdistribu.coords["time"].values[itime_val]}')
        axs[ix2].set_xlabel(f"${f_spec.coords['v_x'].name}$", fontsize=12)
        axs[ix2].set_ylabel('fdistribu', fontsize=12)
        axs[ix2].legend()

        # plot 3: f(x) for 3 v_x values
        ix3 = 3*ispec+2
        for iivx, ivx in enumerate(ivx_array):
            axs[ix3].semilogy(f_spec.coords['x'],
                              np.abs(f_spec.isel(v_x=ivx)),
                              line_type[iivx],
                              label=f'v_x= {fdistribu.coords["v_x"][ivx].item()}')
        axs[ix3].set_title(f'{f_spec_name} at t = {time_diag}')
        axs[ix3].set_xlabel(f"${f_spec.coords['x'].name}$", fontsize=12)
        axs[ix3].set_ylabel('fdistribu', fontsize=12)
        axs[ix3].legend()

    if output_filename is None:
        output_filename = Path(f'fdistribu_t{time_diag}.png')
    print("Save figure in " / output_filename)
    plt.savefig(output_filename)
