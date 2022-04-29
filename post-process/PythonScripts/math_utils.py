# SPDX-License-Identifier: MIT

'''
General useful mathematic functions
'''

# pylint: disable=import-error

import numpy as np
import xarray as xr


def fft1d(data: xr.DataArray):
    """ 1D fft suited to xarray with regular mesh"""

    assert(data.ndim == 1)
    in_dim = data.dims[0]

    in_coords = data.coords[in_dim].values
    coord_steps = in_coords[1:] - in_coords[:-1]
    assert(abs(max(coord_steps) - min(coord_steps)) < .00001)
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
