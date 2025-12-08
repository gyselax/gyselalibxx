// SPDX-License-Identifier: MIT

#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "splitvlasovsolver.hpp"
#include "vector_field_common.hpp"

SplitVlasovSolver::SplitVlasovSolver(
        IAdvectionSpatial<GeometryVxVyXY, GridX> const& advec_x,
        IAdvectionSpatial<GeometryVxVyXY, GridY> const& advec_y,
        IAdvectionVelocity<GeometryVxVyXY, GridVx> const& advec_vx,
        IAdvectionVelocity<GeometryVxVyXY, GridVy> const& advec_vy)
    : m_advec_x(advec_x)
    , m_advec_y(advec_y)
    , m_advec_vx(advec_vx)
    , m_advec_vy(advec_vy)
{
}

DFieldSpVxVyXY SplitVlasovSolver::operator()(
        DFieldSpVxVyXY const allfdistribu,
        DVectorConstFieldXY const electric_field,
        double const dt) const
{
    m_advec_x(allfdistribu, dt / 2);
    m_advec_y(allfdistribu, dt / 2);
    m_advec_vx(allfdistribu, ddcHelper::get<X>(electric_field), dt / 2);
    m_advec_vy(allfdistribu, ddcHelper::get<Y>(electric_field), dt);
    m_advec_vx(allfdistribu, ddcHelper::get<X>(electric_field), dt / 2);
    m_advec_y(allfdistribu, dt / 2);
    m_advec_x(allfdistribu, dt / 2);

    return allfdistribu;
}
