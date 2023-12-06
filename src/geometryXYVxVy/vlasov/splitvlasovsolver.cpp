// SPDX-License-Identifier: MIT

#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "splitvlasovsolver.hpp"

SplitVlasovSolver::SplitVlasovSolver(
        IAdvectionSpatial<GeometryXYVxVy, IDimX> const& advec_x,
        IAdvectionSpatial<GeometryXYVxVy, IDimY> const& advec_y,
        IAdvectionVelocity<GeometryXYVxVy, IDimVx> const& advec_vx,
        IAdvectionVelocity<GeometryXYVxVy, IDimVy> const& advec_vy)
    : m_advec_x(advec_x)
    , m_advec_y(advec_y)
    , m_advec_vx(advec_vx)
    , m_advec_vy(advec_vy)
{
}

DSpanSpXYVxVy SplitVlasovSolver::operator()(
        DSpanSpXYVxVy const allfdistribu_host,
        DViewXY const electric_field_x_host,
        DViewXY const electric_field_y_host,
        double const dt) const
{
    auto allfdistribu_alloc
            = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), allfdistribu_host);
    ddc::ChunkSpan allfdistribu = allfdistribu_alloc.span_view();
    auto electric_field_x_alloc = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), electric_field_x_host);
    ddc::ChunkSpan electric_field_x = electric_field_x_alloc.span_view();
    auto electric_field_y_alloc = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), electric_field_y_host);
    ddc::ChunkSpan electric_field_y = electric_field_y_alloc.span_view();

    m_advec_x(allfdistribu, dt / 2);
    m_advec_y(allfdistribu, dt / 2);
    m_advec_vx(allfdistribu, electric_field_x, dt / 2);
    m_advec_vy(allfdistribu, electric_field_y, dt);
    m_advec_vx(allfdistribu, electric_field_x, dt / 2);
    m_advec_y(allfdistribu, dt / 2);
    m_advec_x(allfdistribu, dt / 2);

    ddc::deepcopy(allfdistribu_host, allfdistribu);
    return allfdistribu_host;
}
