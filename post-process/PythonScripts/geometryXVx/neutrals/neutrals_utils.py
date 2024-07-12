"""
SPDX-License-Identifier: MIT
Utility functions for geometryMX post-process scripts. 
"""

import numpy as np
import xarray as xr

def get_charge_exchange_rate(ds, T_i):
    """
    Returns the charge exchange (CX) rate computed from the rates or reading the constant value from PDI.
    """
    try:
        if ds['charge_exchange_coefficients'].size > 0:
            T_i_log10 = np.log10(T_i)
            cx_coefficients = ds['charge_exchange_coefficients']
            cx_rate_log10 = np.polyval(cx_coefficients, T_i_log10)
            norm_coeff_rate = ds['norm_coeff_rate_neutrals'].values
            k_cx = 10**cx_rate_log10 * norm_coeff_rate
            return k_cx
    except KeyError as e:
        raise KeyError(f"{e} is missing from the dataset.") from e

    try:
        return ds['charge_exchange_constant_rate']
    except KeyError:
        raise KeyError("charge_exchange_constant_rate is missing from the dataset.") from e

def get_ionization_rate(ds, T_e, n_e):
    """
    Returns the ionization rate computed from the rates or reading the constant value from PDI.
    """
    try:
        if ds['ionization_slope_coefficients'].size > 0:
            T_e_log10 = np.log10(T_e)
            density_e_log10 = np.log10(n_e)
            i_slope_coefficient = ds['ionization_slope_coefficients']
            i_intercept_coefficient = ds['ionization_intercept_coefficients']
            norm_coeff_rate = ds['norm_coeff_rate_neutrals'].values

            coefficients_size_proxy = len(i_slope_coefficient) - 1

            i_rate_log10 = xr.zeros_like(density_e_log10, dtype=float)

            for i in range(coefficients_size_proxy):
                polynomial_cs = i_slope_coefficient[i] * density_e_log10 + i_intercept_coefficient[i]
                i_rate_log10 += polynomial_cs * np.power(T_e_log10, i)

            k_i = np.power(10, i_rate_log10) * norm_coeff_rate
            return k_i
    except KeyError as e:
        raise KeyError(f"{e} is missing from the dataset.") from e

    try:
        return ds['ionization_constant_rate']
    except KeyError:
        raise KeyError("charge_exchange_constant_rate is missing from the dataset.") from e

def get_recombination_rate(ds, T_e, n_e):
    """
    Returns the recombination rate computed from the rates or reading the constant value from PDI.
    """
    try:
        if ds['recombination_slope_coefficients'].size > 0:
            T_e_log10 = np.log10(T_e)
            density_e_log10 = np.log10(n_e)
            norm_coeff_rate = ds['norm_coeff_rate_neutrals'].values
            r_slope_coefficient = ds['recombination_slope_coefficients']
            r_intercept_coefficient = ds['recombination_intercept_coefficients']

            coefficients_size_proxy = len(r_slope_coefficient) - 1

            r_rate_log10 = xr.zeros_like(density_e_log10, dtype=float)

            for i in range(coefficients_size_proxy):
                polynomial_cs = r_slope_coefficient[i] * density_e_log10 + r_intercept_coefficient[i]
                r_rate_log10 += polynomial_cs * np.power(T_e_log10, i)

            k_r = np.power(10, r_rate_log10) * norm_coeff_rate
            return k_r
    except KeyError as e:
        raise KeyError(f"{e} is missing from the dataset.") from e

    try:
        return ds['recombination_constant_rate']
    except KeyError:
        raise KeyError("charge_exchange_constant_rate is missing from the dataset.") from e
