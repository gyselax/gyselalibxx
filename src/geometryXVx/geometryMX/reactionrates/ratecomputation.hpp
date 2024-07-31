// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_aliases.hpp"
#include "ireactionrate.hpp"

namespace detail {

template <typename CoeffView>
DFieldSpX compute_rate_electron_density_temperature(
        Kokkos::DefaultExecutionSpace exec_space,
        DFieldSpX rate,
        DConstFieldSpX density,
        DConstFieldSpX temperature,
        CoeffView slope_coefficients,
        CoeffView intercept_coefficients,
        double norm_coeff_rate)
{
    ddc::parallel_for_each(
            exec_space,
            get_idx_range(rate),
            KOKKOS_LAMBDA(IdxSpX const ifspx) {
                double temperature_log10
                        = Kokkos::log10(temperature(ielec(), ddc::select<GridX>(ifspx)));
                double density_log10 = Kokkos::log10(density(ielec(), ddc::select<GridX>(ifspx)));

                double rate_log10 = 0.0;
                for (int i = 0; i < slope_coefficients.extent(0); ++i) {
                    double polynomial_cs
                            = slope_coefficients(i) * density_log10 + intercept_coefficients(i);
                    rate_log10 += polynomial_cs * Kokkos::pow(temperature_log10, i);
                }
                rate(ifspx) = Kokkos::pow(10, rate_log10) * norm_coeff_rate;
            });
    return rate;
}

template <typename CoeffView>
DFieldSpX compute_rate_ion_temperature(
        Kokkos::DefaultExecutionSpace exec_space,
        DFieldSpX rate,
        DConstFieldSpX temperature,
        CoeffView polynomial_coefficients,
        double norm_coeff_rate)
{
    ddc::parallel_for_each(
            exec_space,
            get_idx_range(rate),
            KOKKOS_LAMBDA(IdxSpX const ifspx) {
                double temperature_log10
                        = Kokkos::log10(temperature(ielec(), ddc::select<GridX>(ifspx)));

                double rate_log10 = 0.0;
                for (int i = 0; i < polynomial_coefficients.extent(0); ++i) {
                    rate_log10 += polynomial_coefficients[i] * Kokkos::pow(temperature_log10, i);
                }
                rate(ifspx) = Kokkos::pow(10, rate_log10) * norm_coeff_rate;
            });
    return rate;
}
} // namespace detail
