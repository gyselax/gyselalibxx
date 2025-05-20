// SPDX-License-Identifier: MIT

#include "mpisplithybridvlasovsolver.hpp"

MpiSplitHybridVlasovSolver::MpiSplitHybridVlasovSolver(
        IAdvectionSpatial<GeometryVxVyXY, GridX> const& advec_x,
        IAdvectionSpatial<GeometryVxVyXY, GridY> const& advec_y,
        IAdvectionVelocity<GeometryXYVxVy, GridVx> const& advec_vx,
        IAdvectionVelocity<GeometryXYVxVy, GridVy> const& advec_vy,
        IAdvectionVelocityRot2D<GeometryXYVxVy, GridVx> const& advec_2d_rot_vx,
        IAdvectionVelocityRot2D<GeometryXYVxVy, GridVy> const& advec_2d_rot_vy,
        MPITransposeAllToAll<X2DSplit, V2DSplit> const& transpose)
    : m_advec_x(advec_x)
    , m_advec_y(advec_y)
    , m_advec_vx(advec_vx)
    , m_advec_vy(advec_vy)
    , m_advec_2d_rot_vx(advec_2d_rot_vx)
    , m_advec_2d_rot_vy(advec_2d_rot_vy)
    , m_transpose(transpose)
{
}

DFieldSpVxVyXY MpiSplitHybridVlasovSolver::operator()(
        DFieldSpVxVyXY const allfdistribu_v2Dsplit,
        DFieldXY magnetic_field_z,
        DFieldXY frame_shift_x,
        DFieldXY frame_shift_y,
        double const dt) const
{
    IdxRangeSpVxVyXY idxrange_v2Dsplit(m_transpose.get_local_idx_range<V2DSplit>());
    IdxRangeSpXYVxVy idxrange_x2Dsplit(m_transpose.get_local_idx_range<X2DSplit>());
    IdxRangeXY idx_range_xy_v2Dsplit(idxrange_v2Dsplit);
    DFieldMemSpXYVxVy allfdistribu_x2Dsplit_alloc(idxrange_x2Dsplit);
    DFieldSpXYVxVy allfdistribu_x2Dsplit = get_field(allfdistribu_x2Dsplit_alloc);

    // Swap to vxvy contiguous layout
    m_transpose(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_x2Dsplit,
            get_const_field(allfdistribu_v2Dsplit));
    // Advect in velocity dimensions: rotation by exact splittings
    m_advec_2d_rot_vx(allfdistribu_x2Dsplit, get_const_field(magnetic_field_z), get_const_field(frame_shift_x), get_const_field(frame_shift_y), dt);
    m_advec_2d_rot_vy(allfdistribu_x2Dsplit, get_const_field(magnetic_field_z), get_const_field(frame_shift_x), get_const_field(frame_shift_y), dt);
    m_advec_2d_rot_vx(allfdistribu_x2Dsplit, get_const_field(magnetic_field_z), get_const_field(frame_shift_x), get_const_field(frame_shift_y), dt);
    
    // Swap to xy contiguous layout
    m_transpose(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_v2Dsplit,
            get_const_field(allfdistribu_x2Dsplit));

    return allfdistribu_v2Dsplit;
}


DFieldSpVxVyXY MpiSplitHybridVlasovSolver::operator()(
        DFieldSpVxVyXY const allfdistribu_v2Dsplit,
        double const dt) const
{
    // Advect in spatial dimensions
    m_advec_x(allfdistribu_v2Dsplit, dt);
    m_advec_y(allfdistribu_v2Dsplit, dt);

    return allfdistribu_v2Dsplit;
}
