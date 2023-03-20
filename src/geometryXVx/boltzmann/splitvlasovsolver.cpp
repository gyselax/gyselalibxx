// SPDX-License-Identifier: MIT

#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "splitvlasovsolver.hpp"

SplitVlasovSolver::SplitVlasovSolver(IAdvectionX const& advec_x, IAdvectionVx const& advec_vx)
    : m_advec_x(advec_x)
    , m_advec_vx(advec_vx)
{
}


DSpanSpXVx SplitVlasovSolver::operator()(
        DSpanSpXVx const allfdistribu,
        DViewX const electric_field,
        double const dt) const
{
    m_advec_x(allfdistribu, dt / 2);
    m_advec_vx(allfdistribu, electric_field, dt);
    m_advec_x(allfdistribu, dt / 2);
    return allfdistribu;
}
