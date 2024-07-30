// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "constantfluidinitialization.hpp"

ConstantFluidInitialization::ConstantFluidInitialization(host_t<DViewSpM> moments)
    : m_moments_alloc(moments.domain())
{
    ddc::parallel_deepcopy(m_moments_alloc, moments);
}

DSpanSpMX ConstantFluidInitialization::operator()(DSpanSpMX const fluid_moments) const
{
    ddc::ChunkSpan moments(m_moments_alloc.span_view());
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            fluid_moments.domain(),
            KOKKOS_LAMBDA(IndexSpMX const ispmx) {
                IndexSpM ispm(ddc::select<Species, IDimM>(ispmx));
                fluid_moments(ispmx) = moments(ispm);
            });
    return fluid_moments;
}
