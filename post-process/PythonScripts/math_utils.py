# SPDX-License-Identifier: MIT

'''
General useful mathematic functions
'''

# pylint: disable=import-error

import numpy as np
import xarray as xr


def fft1d(data: xr.DataArray):
    """ 1D fft suited to xarray with regular mesh"""

    assert data.ndim == 1
    in_dim = data.dims[0]

    in_coords = data.coords[in_dim].values
    coord_steps = in_coords[1:] - in_coords[:-1]
    assert abs(max(coord_steps) - min(coord_steps)) < .00001
    coord_step = coord_steps[0]

    result = np.fft.fft(data)
    freqs = (2.*np.pi/coord_step) * np.fft.fftfreq(data.size)
    result = np.fft.fftshift(result)
    freqs = np.fft.fftshift(freqs)

    out_dim = in_dim+"_freq"
    out_coords = dict(data.coords)
    del out_coords[in_dim]
    out_coords[out_dim] = freqs
    result = xr.DataArray(result, dims=[out_dim], coords=out_coords)

    return result


def differentiate(data, coord_name, **opt_args):
    """
    differentiates data along the coord_name coordinate

    Parameters
    ----------
    data: xarray.DataArray
        data to integrate

    coord_name: coordinate along which to differentiate
    """

    if coord_name not in data.coords:
        raise ValueError(f"Coordinate {coord_name} is not a coordinate of data, {data.coords}")

    return data.differentiate(coord_name, edge_order=2)


def compute_moment(data, coord_name, moment_order, coeff=1, **opt_args):
    """
    Integrates coeff*data*coord^moment_order along the coord_name coordinate

    Parameters
    ----------
    data: xarray.DataArray
        data to integrate

    coord_name: coordinate along which to integrate
    """

    if coord_name not in data.coords:
        raise ValueError(f"Coordinate {coord_name} is not a coordinate of data, {data.coords}")

    if moment_order == 0:
        return data.integrate(coord_name)

    coordinate = data.coords[coord_name]
    data_moment = coeff*data*coordinate.values**moment_order
    return data_moment.integrate(coord_name)
