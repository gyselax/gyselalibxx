// SPDX-License-Identifier: MIT

#include "nullfluidsolver.hpp"

NullFluidSolver::NullFluidSolver(IdxRangeSp const& dom_fluidsp)
{
    // charged fluid species is not allowed for now
    for (IdxSp const isp : dom_fluidsp) {
        assert(charge(isp) == 0.);
    }
}

DSpanSpMX NullFluidSolver::operator()(
        DSpanSpMX const fluid_moments,
        DViewSpXVx const allfdistribu,
        DViewX const efield,
        double const dt) const
{
    return fluid_moments;
}
