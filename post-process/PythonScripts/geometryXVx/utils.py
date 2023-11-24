# SPDX-License-Identifier: MIT

import os
import yaml


def get_simulation_parameters(path, filename):
    '''Parses a yaml file containing the simulation parameters
    '''
    with open(os.path.join(path, filename)) as file:
        data = yaml.safe_load(file)

        if len(data['SpeciesInfo']) != 2:
            raise ValueError('The number of species is not two')

        if data['SpeciesInfo'][0]['charge'] * data['SpeciesInfo'][1]['charge'] > 0:
            raise ValueError('The two species do not have opposite charges')

    return data


def compute_fluid_velocity(density, particle_flux):
    return particle_flux / density


def compute_pressure(density, particle_flux, momentum_flux):
    fluid_velocity = particle_flux / density
    return momentum_flux - particle_flux*fluid_velocity


def compute_temperature(density, particle_flux, momentum_flux):
    return compute_pressure(density, particle_flux, momentum_flux) / density


def compute_krook_sink_constant(diskstore):
    nu = diskstore['krook_sink_constant_amplitude']
    mask = diskstore['krook_sink_constant_mask']
    ftarget = diskstore['krook_sink_constant_ftarget']
    return -nu*mask*(diskstore['fdistribu']-ftarget)


def compute_kinetic_source(diskstore):
    return diskstore['kinetic_source_amplitude'] \
            * diskstore['kinetic_source_spatial_extent'] * diskstore['kinetic_source_velocity_shape']
