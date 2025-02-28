// SPDX-License-Identifier: MIT

#include "mpisplitvlasovsolver.hpp"

MpiSplitVlasovSolver::MpiSplitVlasovSolver(
        IAdvectionSpatial<GeometryXYVxVy, GridX> const& advec_x,
        IAdvectionSpatial<GeometryXYVxVy, GridY> const& advec_y,
        IAdvectionVelocity<GeometryXYVxVy, GridVx> const& advec_vx,
        IAdvectionVelocity<GeometryXYVxVy, GridVy> const& advec_vy,
        MPITransposeAllToAll<X2DSplit, V2DSplit> const& transpose)
    : m_advec_x(advec_x)
    , m_advec_y(advec_y)
    , m_advec_vx(advec_vx)
    , m_advec_vy(advec_vy)
    , m_transpose(transpose)
{
}

DFieldSpXYVxVy MpiSplitVlasovSolver::operator()(
        DFieldSpXYVxVy const allfdistribu_x2Dsplit,
        DConstFieldXY const electric_field_x,
        DConstFieldXY const electric_field_y,
        double const dt) const
{
    IdxRangeSpVxVyXY idxrange_v2Dsplit(transpose.get_local_idx_range<V2DSplit>());
    IdxRangeXY idx_range_xy_v2Dsplit(idxrange_v2Dsplit);
    DFieldMemSpVxVyXY const allfdistribu_v2Dsplit_alloc(idxrange_v2Dsplit);
    DFieldSpVxVyXY const allfdistribu_v2Dsplit(allfdistribu_v2Dsplit_alloc);
    m_advec_vx(allfdistribu_x2Dsplit, electric_field_x[idx_range_xy_v2Dsplit], dt / 2);
    m_advec_vy(allfdistribu_x2Dsplit, electric_field_y[idx_range_xy_v2Dsplit], dt / 2);
    transpose(allfdistribu_v2Dsplit, allfdistribu_x2Dsplit);
    m_advec_x(allfdistribu_v2Dsplit, dt / 2);
    m_advec_y(allfdistribu_v2Dsplit, dt);
    m_advec_x(allfdistribu_v2Dsplit, dt / 2);
    transpose(allfdistribu_x2Dsplit, allfdistribu_v2Dsplit);
    m_advec_vx(allfdistribu_x2Dsplit, electric_field_x[idx_range_xy_v2Dsplit], dt / 2);
    m_advec_vy(allfdistribu_x2Dsplit, electric_field_y[idx_range_xy_v2Dsplit], dt / 2);

    return allfdistribu_x2Dsplit;
}
