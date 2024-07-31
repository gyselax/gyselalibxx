// SPDX-License-Identifier: MIT

#include "ionization.hpp"
#include "ratecomputation.hpp"

IonizationRate::IonizationRate(double const norm_coeff_rate)
    : m_slope_coefficient("slope_coefficient")
    , m_intercept_coefficient("intercept_coefficient")
    , m_norm_coeff_rate(norm_coeff_rate)
{
    double slope_coefficient_alloc[s_coefficients_size]
            = {0.06478416786452716,
               -0.06283835812957501,
               0.0533047026400489,
               -0.02412714188219726,
               0.005377140897108906,
               -0.00047180128893166466};
    double intercept_coefficient_alloc[s_coefficients_size]
            = {0.053216546399122204,
               1.4939944926753081,
               -1.8156125022306027,
               1.269858057577872,
               -0.468107496181505,
               0.06406763661421713};

    DKokkosView_h<s_coefficients_size> slope_coefficient_host(slope_coefficient_alloc);
    DKokkosView_h<s_coefficients_size> intercept_coefficient_host(intercept_coefficient_alloc);
    Kokkos::deep_copy(m_slope_coefficient, slope_coefficient_host);
    Kokkos::deep_copy(m_intercept_coefficient, intercept_coefficient_host);

    PDI_multi_expose(
            "i_rate_coeff_pol_expose",
            "ionization_slope_coefficients",
            slope_coefficient_host.data(),
            PDI_OUT,
            "ionization_intercept_coefficients",
            intercept_coefficient_host.data(),
            PDI_OUT,
            NULL);
}

DFieldSpX IonizationRate::operator()(
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
