// SPDX-License-Identifier: MIT

#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "splitvlasovsolver.hpp"

SplitVlasovSolver::SplitVlasovSolver(
        IAdvectionSpatial<GeometryXVx, IDimX> const& advec_x,
        IAdvectionVelocity<GeometryXVx, IDimVx> const& advec_vx)
    : m_advec_x(advec_x)
    , m_advec_vx(advec_vx)
{
}

device_t<DSpanSpXVx> SplitVlasovSolver::operator()(
        device_t<DSpanSpXVx> const allfdistribu,
        device_t<DViewX> const electric_field,
        double const dt) const
{
    m_advec_x(allfdistribu, dt / 2);
    m_advec_vx(allfdistribu, electric_field, dt);
    m_advec_x(allfdistribu, dt / 2);
    return allfdistribu;
}
