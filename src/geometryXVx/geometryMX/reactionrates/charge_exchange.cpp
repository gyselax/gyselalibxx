// SPDX-License-Identifier: MIT

#include <pdi.h>

#include "charge_exchange.hpp"
#include "ratecomputation.hpp"
#include "species_info.hpp"

ChargeExchangeRate::ChargeExchangeRate(double const norm_coeff_rate)
    : m_cx_coefficients("chargexchange_coefficients")
    , m_norm_coeff_rate(norm_coeff_rate)
{
    double cx_coefficients_alloc[s_coefficients_size]
            = {0.25125030648153424,
               0.37349850950319474,
               -0.009966677207985306,
               0.017922185310747015,
               -0.010699854255929037};

    DKokkosView_h<s_coefficients_size> cx_coefficients_host(cx_coefficients_alloc);
    Kokkos::deep_copy(m_cx_coefficients, cx_coefficients_host);
    PDI_multi_expose(
            "cx_rate_coeff_pol_expose",
            "charge_exchange_coefficients",
            cx_coefficients_host.data(),
            PDI_OUT,
            NULL);
}

IdxSp find_ion(IdxRangeSp const idx_range_kinsp)
{
    bool ion_found = false;
    IdxSp iion;
    for (IdxSp const isp : idx_range_kinsp) {
        if (charge(isp) > 0.) {
            ion_found = true;
            iion = isp;
        }
    }
    if (!ion_found) {
        throw std::runtime_error("ion not found");
    }
    assert(idx_range_kinsp.size() == 2);

    return iion;
}

DFieldSpX ChargeExchangeRate::operator()(
        DFieldSpX rate,
        DConstFieldSpX density,
        DConstFieldSpX temperature) const
{
    return detail::compute_rate_ion_temperature(
            Kokkos::DefaultExecutionSpace(),
            rate,
            temperature,
            m_cx_coefficients,
            m_norm_coeff_rate);
}
