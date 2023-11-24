"""
SPDX-License-Identifier: MIT
Utility functions for geometryXVx post-process scripts. 
"""

import os
import yaml
import xarray as xr


def get_simulation_parameters(path, filename):
    '''Parses a yaml file containing the simulation parameters
    '''
    with open(os.path.join(path, filename), encoding='utf8') as file:
        data = yaml.safe_load(file)

        if len(data['SpeciesInfo']) != 2:
            raise ValueError('The number of species is not two')

        if data['SpeciesInfo'][0]['charge'] * data['SpeciesInfo'][1]['charge'] > 0:
            raise ValueError('The two species do not have opposite charges')

    return data


def compute_fluid_velocity(density, particle_flux):
    '''Computes the fluid velocity moment
    '''
    return particle_flux / density


def compute_pressure(density, particle_flux, momentum_flux):
    '''Computes the pressure moment
    '''
    fluid_velocity = particle_flux / density
    return momentum_flux - particle_flux*fluid_velocity


def compute_temperature(density, particle_flux, momentum_flux):
    '''Computes the temperature moment
    '''
    return compute_pressure(density, particle_flux, momentum_flux) / density


def compute_krook_sink_constant(diskstore):
    '''Computes the constant krook sink expression
    '''
    try:
        nu = diskstore['krook_sink_constant_amplitude']
        mask = diskstore['krook_sink_constant_mask']
        ftarget = diskstore['krook_sink_constant_ftarget']

    except KeyError as e:
        print('Info: no constant krook sink in simulation:', e)
        return xr.zeros_like(diskstore['fdistribu'])

    return -nu*mask*(diskstore['fdistribu']-ftarget)


def compute_krook_sink_adaptive(diskstore, density):
    '''Computes the adaptive krook sink expression
    '''
    try:
        mask = diskstore['krook_sink_adaptive_mask']
        ftarget = diskstore['krook_sink_adaptive_ftarget']
        nu_ions = diskstore['krook_sink_adaptive_amplitude']
        density_target = diskstore['krook_sink_adaptive_density']
        nu = nu_ions*(density.sel(species='ions')-density_target) \
            / (density-density_target)

    except KeyError as e:
        print('Info: no adaptive krook sink in simulation:', e)
        return xr.zeros_like(diskstore['fdistribu'])

    return -nu*mask*(diskstore['fdistribu']-ftarget)


def compute_kinetic_source(diskstore):
    '''Computes the kinetic source expression
    '''
    try:
        return diskstore['kinetic_source_amplitude'] \
                * diskstore['kinetic_source_spatial_extent'] \
                * diskstore['kinetic_source_velocity_shape']
    except KeyError as e:
        print('Info: no kinetic source in simulation:', e)
        return xr.zeros_like(diskstore['fdistribu'])

