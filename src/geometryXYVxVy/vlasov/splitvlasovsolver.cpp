// SPDX-License-Identifier: MIT

#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "splitvlasovsolver.hpp"

SplitVlasovSolver::SplitVlasovSolver(
        IAdvectionSpatial<GeometryXYVxVy, GridX> const& advec_x,
        IAdvectionSpatial<GeometryXYVxVy, GridY> const& advec_y,
        IAdvectionVelocity<GeometryXYVxVy, GridVx> const& advec_vx,
        IAdvectionVelocity<GeometryXYVxVy, GridVy> const& advec_vy)
    : m_advec_x(advec_x)
    , m_advec_y(advec_y)
    , m_advec_vx(advec_vx)
    , m_advec_vy(advec_vy)
{
}

DFieldSpXYVxVy SplitVlasovSolver::operator()(
        DFieldSpXYVxVy const allfdistribu,
        DConstFieldXY const electric_field_x,
        DConstFieldXY const electric_field_y,
        double const dt) const
{
    m_advec_x(allfdistribu, dt / 2);
    m_advec_y(allfdistribu, dt / 2);
    m_advec_vx(allfdistribu, electric_field_x, dt / 2);
    m_advec_vy(allfdistribu, electric_field_y, dt);
    m_advec_vx(allfdistribu, electric_field_x, dt / 2);
    m_advec_y(allfdistribu, dt / 2);
    m_advec_x(allfdistribu, dt / 2);

    return allfdistribu;
}
