// SPDX-License-Identifier: MIT

#include "mpisplitvlasovsolver.hpp"

MpiSplitVlasovSolver::MpiSplitVlasovSolver(
        IAdvectionSpatial<GeometryVxVyVzX, GridX> const& advec_x,
        IAdvectionVelocity<GeometryXVxVyVz, GridVx> const& advec_vx,
        IAdvectionVelocity<GeometryXVxVyVz, GridVy> const& advec_vy,
        IAdvectionVelocity<GeometryXVxVyVz, GridVz> const& advec_vz,
        MPITransposeAllToAll<X1DSplit, V3DSplit> const& transpose)
    : m_advec_x(advec_x)
    , m_advec_vx(advec_vx)
    , m_advec_vy(advec_vy)
    , m_advec_vz(advec_vz)
    , m_transpose(transpose)
{
}

DFieldSpVxVyVzX MpiSplitVlasovSolver::operator()(
        DFieldSpVxVyVzX const allfdistribu_v3Dsplit,
        DConstFieldX const electric_field_x,
        DConstFieldX const electric_field_y,
        DConstFieldX const electric_field_z,
        double const dt) const
{
    IdxRangeSpVxVyVzX idxrange_v3Dsplit(m_transpose.get_local_idx_range<V3DSplit>());
    IdxRangeSpXVxVyVz idxrange_x1Dsplit(m_transpose.get_local_idx_range<X1DSplit>());
    IdxRangeX idx_range_x_v3Dsplit(idxrange_v3Dsplit);
    DFieldMemSpXVxVyVz allfdistribu_x1Dsplit_alloc(idxrange_x1Dsplit);
    DFieldSpXVxVyVz allfdistribu_x1Dsplit = get_field(allfdistribu_x1Dsplit_alloc);

    // Create contiguous memory space to contain the relevant section of the electric field
    DFieldMemX local_electric_field_x(idx_range_x_v3Dsplit);
    DFieldMemX local_electric_field_y(idx_range_x_v3Dsplit);
    DFieldMemX local_electric_field_z(idx_range_x_v3Dsplit);
    ddc::parallel_deepcopy(
            get_field(local_electric_field_x),
            electric_field_x[idx_range_x_v3Dsplit]);
    ddc::parallel_deepcopy(
            get_field(local_electric_field_y),
            electric_field_y[idx_range_x_v3Dsplit]);
    ddc::parallel_deepcopy(
            get_field(local_electric_field_z),
            electric_field_z[idx_range_x_v3Dsplit]);

    // Advect in spatial dimensions
    m_advec_x(allfdistribu_v3Dsplit, dt / 2);
    // Swap to vxvyvz contiguous layout
    m_transpose(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_x1Dsplit,
            get_const_field(allfdistribu_v3Dsplit));
    // Advect in velocity dimensions
    m_advec_vx(allfdistribu_x1Dsplit, get_const_field(local_electric_field_x), dt / 2);
    m_advec_vy(allfdistribu_x1Dsplit, get_const_field(local_electric_field_y), dt / 2);
    m_advec_vz(allfdistribu_x1Dsplit, get_const_field(local_electric_field_z), dt);
    m_advec_vy(allfdistribu_x1Dsplit, get_const_field(local_electric_field_y), dt / 2);
    m_advec_vx(allfdistribu_x1Dsplit, get_const_field(local_electric_field_x), dt / 2);
    // Swap to xy contiguous layout
    m_transpose(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_v3Dsplit,
            get_const_field(allfdistribu_x1Dsplit));
    // Advect in spatial dimensions
    m_advec_x(allfdistribu_v3Dsplit, dt / 2);

    return allfdistribu_v3Dsplit;
}
