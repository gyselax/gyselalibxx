// SPDX-License-Identifier: MIT

#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "splitvlasovsolver.hpp"

SplitVlasovSolver::SplitVlasovSolver(
        IAdvectionSpatial<GeometryVxVyVzX, GridX> const& advec_x,
        IAdvectionVelocity<GeometryVxVyVzX, GridVx> const& advec_vx,
        IAdvectionVelocity<GeometryVxVyVzX, GridVy> const& advec_vy,
        IAdvectionVelocity<GeometryVxVyVzX, GridVz> const& advec_vz)
    : m_advec_x(advec_x)
    , m_advec_vx(advec_vx)
    , m_advec_vy(advec_vy)
    , m_advec_vz(advec_vz)
{
}

DFieldSpVxVyVzX SplitVlasovSolver::operator()(
        DFieldSpVxVyVzX const allfdistribu,
        DConstFieldX const electric_field_x,
        DConstFieldX const electric_field_y,
        DConstFieldX const electric_field_z,
        double const dt) const
{
    m_advec_x(allfdistribu, dt / 2);
    m_advec_vx(allfdistribu, electric_field_x, dt / 2);
    m_advec_vy(allfdistribu, electric_field_y, dt / 2);
    m_advec_vz(allfdistribu, electric_field_z, dt);
    m_advec_vy(allfdistribu, electric_field_y, dt / 2);
    m_advec_vx(allfdistribu, electric_field_x, dt / 2);
    m_advec_x(allfdistribu, dt / 2);

    return allfdistribu;
}
