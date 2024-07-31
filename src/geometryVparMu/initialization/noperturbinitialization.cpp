// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "noperturbinitialization.hpp"


NoPerturbInitialization::NoPerturbInitialization(DConstFieldSpVparMu fequilibrium)
    : m_fequilibrium(fequilibrium)
{
}

DFieldSpVparMu NoPerturbInitialization::operator()(DFieldSpVparMu const allfdistribu) const
{
    IdxRangeSpVparMu const idxrange_spvparmu = get_idx_range(allfdistribu);

    DConstFieldSpVparMu fequilibrium_proxy = get_field(m_fequilibrium);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idxrange_spvparmu,
            KOKKOS_LAMBDA(IdxSpVparMu const ispvparmu) {
                double fdistribu_val = fequilibrium_proxy(ispvparmu);
                if (fdistribu_val < 1.e-60) {
                    fdistribu_val = 1.e-60;
                }
                allfdistribu(ispvparmu) = fdistribu_val;
            });
    return allfdistribu;
}
