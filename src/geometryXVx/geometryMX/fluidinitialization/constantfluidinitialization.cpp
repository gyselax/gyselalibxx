// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "constantfluidinitialization.hpp"

ConstantFluidInitialization::ConstantFluidInitialization(host_t<DConstFieldSpMom> moments)
    : m_moments_alloc(get_idx_range(moments))
{
    ddc::parallel_deepcopy(m_moments_alloc, moments);
}

DFieldSpMomX ConstantFluidInitialization::operator()(DFieldSpMomX const fluid_moments) const
{
    DConstFieldSpMom moments(get_field(m_moments_alloc));
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(fluid_moments),
            KOKKOS_LAMBDA(IdxSpMomX const ispmx) {
                IdxSpMom ispm(ddc::select<Species, GridMom>(ispmx));
                fluid_moments(ispmx) = moments(ispm);
            });
    return fluid_moments;
}
