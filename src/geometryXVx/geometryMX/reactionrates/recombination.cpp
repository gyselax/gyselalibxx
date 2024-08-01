// SPDX-License-Identifier: MIT

#include <species_info.hpp>

#include "ratecomputation.hpp"
#include "recombination.hpp"

RecombinationRate::RecombinationRate(double const norm_coeff_rate)
    : m_slope_coefficient("slope_coefficient")
    , m_intercept_coefficient("intercept_coefficient")
    , m_norm_coeff_rate(norm_coeff_rate)
{
    double slope_coefficient_alloc[s_coefficients_size]
            = {0.03849324670136343, -0.03713312941661376};
    double intercept_coefficient_alloc[s_coefficients_size]
            = {-5.117013348529228, -1.1497231886316353};

    DKokkosView_h<s_coefficients_size> slope_coefficient_host(slope_coefficient_alloc);
    DKokkosView_h<s_coefficients_size> intercept_coefficient_host(intercept_coefficient_alloc);
    Kokkos::deep_copy(m_slope_coefficient, slope_coefficient_host);
    Kokkos::deep_copy(m_intercept_coefficient, intercept_coefficient_host);

    PDI_multi_expose(
            "r_rate_coeff_pol_expose",
            "recombination_slope_coefficients",
            slope_coefficient_host.data(),
            PDI_OUT,
            "recombination_intercept_coefficients",
            intercept_coefficient_host.data(),
            PDI_OUT,
            NULL);
}

DFieldSpX RecombinationRate::operator()(
        DFieldSpX rate,
        DConstFieldSpX density,
        DConstFieldSpX temperature) const
{
    return detail::compute_rate_electron_density_temperature(
            Kokkos::DefaultExecutionSpace(),
            rate,
            density,
            temperature,
            m_slope_coefficient,
            m_intercept_coefficient,
            m_norm_coeff_rate);
}
