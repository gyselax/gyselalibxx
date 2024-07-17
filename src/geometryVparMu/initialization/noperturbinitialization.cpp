// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "noperturbinitialization.hpp"


NoPerturbInitialization::NoPerturbInitialization(DViewSpVparMu fequilibrium)
    : m_fequilibrium(fequilibrium)
{
}

DSpanSpVparMu NoPerturbInitialization::operator()(DSpanSpVparMu const allfdistribu) const
{
    IdxRangeSpVparMu const idxrange_spvparmu = allfdistribu.domain();

    ddc::ChunkSpan fequilibrium_proxy = m_fequilibrium.span_view();
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
