// SPDX-License-Identifier: MIT

#include "mpisplithybridvlasovsolver.hpp"

MpiSplitHybridVlasovSolver::MpiSplitHybridVlasovSolver(
        IAdvectionSpatial<GeometryVxVyVzX, GridX> const& advec_x,
        IAdvectionVelocity<GeometryXVxVyVz, GridVx> const& advec_vx,
        IAdvectionVelocity<GeometryXVxVyVz, GridVy> const& advec_vy,
        IAdvectionVelocity<GeometryXVxVyVz, GridVz> const& advec_vz,
        IAdvectionVelocityRot3D<GeometryXVxVyVz, GridVx, GridVy, GridVz> const& advec_3d_rot,
        MPITransposeAllToAll<X1DSplit, V3DSplit> const& transpose)
    : m_advec_x(advec_x)
    , m_advec_vx(advec_vx)
    , m_advec_vy(advec_vy)
    , m_advec_vz(advec_vz)
    , m_advec_3d_rot(advec_3d_rot)
    , m_transpose(transpose)
{
}

DFieldSpVxVyVzX MpiSplitHybridVlasovSolver::operator()(
        DFieldSpVxVyVzX const allfdistribu_v3Dsplit,
        DFieldX frame_shift_x,
        DFieldX frame_shift_y,
        DFieldX frame_shift_z,
        DFieldX p_parallel_x,
        DFieldX p_parallel_y,
        DFieldX p_parallel_z,
        DFieldX magnetic_x, 
        DFieldX magnetic_y, 
        DFieldX magnetic_z,
        double const dt) const
{
    IdxRangeSpVxVyVzX idxrange_v3Dsplit(m_transpose.get_local_idx_range<V3DSplit>());
    IdxRangeSpXVxVyVz idxrange_x1Dsplit(m_transpose.get_local_idx_range<X1DSplit>());
    IdxRangeX idx_range_x_v3Dsplit(idxrange_v3Dsplit);
    DFieldMemSpXVxVyVz allfdistribu_x1Dsplit_alloc(idxrange_x1Dsplit);
    DFieldSpXVxVyVz allfdistribu_x1Dsplit = get_field(allfdistribu_x1Dsplit_alloc);

    // Create contiguous memory space to contain the relevant section of the electric field
    DFieldMemX local_p_parallel_x(idx_range_x_v3Dsplit);
    DFieldMemX local_p_parallel_y(idx_range_x_v3Dsplit);
    DFieldMemX local_p_parallel_z(idx_range_x_v3Dsplit);

    DFieldMemX local_frame_shift_x(idx_range_x_v3Dsplit);
    DFieldMemX local_frame_shift_y(idx_range_x_v3Dsplit);
    DFieldMemX local_frame_shift_z(idx_range_x_v3Dsplit);

    DFieldMemX local_magnetic_x(idx_range_x_v3Dsplit);
    DFieldMemX local_magnetic_y(idx_range_x_v3Dsplit);
    DFieldMemX local_magnetic_z(idx_range_x_v3Dsplit);

    ddc::parallel_deepcopy(
            get_field(local_p_parallel_x),
            p_parallel_x[idx_range_x_v3Dsplit]);
    ddc::parallel_deepcopy(
            get_field(local_p_parallel_y),
            p_parallel_y[idx_range_x_v3Dsplit]);
    ddc::parallel_deepcopy(
            get_field(local_p_parallel_z),
            p_parallel_z[idx_range_x_v3Dsplit]);

    ddc::parallel_deepcopy(
            get_field(local_frame_shift_x),
            frame_shift_x[idx_range_x_v3Dsplit]);
    ddc::parallel_deepcopy(
            get_field(local_frame_shift_y),
            frame_shift_y[idx_range_x_v3Dsplit]);
    ddc::parallel_deepcopy(
            get_field(local_frame_shift_z),
            frame_shift_z[idx_range_x_v3Dsplit]);

    ddc::parallel_deepcopy(
            get_field(local_magnetic_x),
            magnetic_x[idx_range_x_v3Dsplit]);
    ddc::parallel_deepcopy(
            get_field(local_magnetic_y),
            magnetic_y[idx_range_x_v3Dsplit]);
    ddc::parallel_deepcopy(
            get_field(local_magnetic_z),
            magnetic_z[idx_range_x_v3Dsplit]);

    // Swap to vxvyvz contiguous layout
    m_transpose(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_x1Dsplit,
            get_const_field(allfdistribu_v3Dsplit));
    
    m_advec_3d_rot(allfdistribu_x1Dsplit, get_const_field(local_magnetic_x), get_const_field(local_magnetic_y), get_const_field(local_magnetic_z), get_const_field(local_frame_shift_x), get_const_field(local_frame_shift_y), get_const_field(local_frame_shift_z), dt);

    m_advec_vx(allfdistribu_x1Dsplit, get_const_field(local_p_parallel_x), dt);
    m_advec_vy(allfdistribu_x1Dsplit, get_const_field(local_p_parallel_y), dt);
    m_advec_vz(allfdistribu_x1Dsplit, get_const_field(local_p_parallel_z), dt);
    
    // Swap to x contiguous layout
    m_transpose(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_v3Dsplit,
            get_const_field(allfdistribu_x1Dsplit));
    //std::cout << "check v55555 translation: " << std::endl;

    return allfdistribu_v3Dsplit;

}


DFieldSpVxVyVzX MpiSplitHybridVlasovSolver::operator()(
        DFieldSpVxVyVzX const allfdistribu_v3Dsplit,
        double const dt) const
{
    // Advect in spatial dimensions
    m_advec_x(allfdistribu_v3Dsplit, dt);

    return allfdistribu_v3Dsplit;
}
