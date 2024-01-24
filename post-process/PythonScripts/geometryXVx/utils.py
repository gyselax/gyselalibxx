"""
SPDX-License-Identifier: MIT
Utility functions for geometryXVx post-process scripts. 
"""

import os
import numpy as np
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


def compute_collinter_collision_frequency(diskstore, density, temperature):
    '''Computes the interspecies collision frequency
    '''

    charges = diskstore['fdistribu_charges']
    masses = diskstore['fdistribu_masses']
    charges_pow4 = charges**4
    collfreq_onespec = diskstore['collinter_nustar0'] * charges_pow4 \
        * np.sqrt(masses.sel(species='electrons')/masses) \
        / diskstore['Lx'] * density / np.power(temperature, 1.5)

    density_flipped = density.copy()
    density_flipped.coords['species'] = density.coords['species'][::-1]
    density_ratio = density/density_flipped

    temperature_flipped = temperature.copy()
    temperature_flipped.coords['species'] = temperature.coords['species'][::-1]
    temperature_ratio = temperature/temperature_flipped

    charges_ratio = charges/charges.values[::-1]
    masses_ratio = masses/masses.values[::-1]

    collfreq_bothspecies = collfreq_onespec * np.sqrt(2.) / (charges_ratio * charges_ratio) \
        * (1. + masses_ratio) / density_ratio \
        / np.power(1. + masses_ratio/temperature_ratio, 1.5)

    return collfreq_bothspecies


def compute_collinter_momentum_exchange(diskstore, density, fluid_velocity, temperature):
    '''Computes the interspecies collision operator expression of the 
    momentum exchange term
    '''
    try:
        collfreq = compute_collinter_collision_frequency(diskstore, density, temperature)
        masses = diskstore['fdistribu_masses']
        masses_ratio = masses/masses.values[::-1]
        fluid_velocity_flipped = fluid_velocity.copy()
        fluid_velocity_flipped.coords['species'] = fluid_velocity_flipped.coords['species'][::-1]
        momentum_term = -collfreq*density * (fluid_velocity - np.sqrt(masses_ratio) \
                                    * fluid_velocity_flipped)
        return momentum_term
    except KeyError as e:
        print('Info: no inter_species collisions in simulation:', e)
        return xr.zeros_like(diskstore['electrostatic_potential'])


def compute_collinter_energy_exchange(diskstore, density, temperature):
    '''Computes the interspecies collision operator expression of the 
    energy exchange term
    '''
    try:
        temperature_flipped = temperature.copy()
        temperature_flipped.coords['species'] = temperature_flipped.coords['species'][::-1]

        masses = diskstore['fdistribu_masses']
        masses_flipped = masses.copy()
        masses_flipped.coords['species'] = masses_flipped.coords['species'][::-1]
        ma_on_ma_mb = masses/(masses + masses_flipped)
        collfreq = compute_collinter_collision_frequency(diskstore, density, temperature)

        return -3*collfreq*density*ma_on_ma_mb*(temperature - temperature_flipped)

    except KeyError as e:
        print('Info: no inter_species collisions in simulation:', e)
        return xr.zeros_like(diskstore['electrostatic_potential'])
