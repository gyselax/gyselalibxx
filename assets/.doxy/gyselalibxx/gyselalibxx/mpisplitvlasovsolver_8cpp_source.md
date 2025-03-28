

# File mpisplitvlasovsolver.cpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**vlasov**](dir_0a9688649b1824bbfb2c211b845ba732.md) **>** [**mpisplitvlasovsolver.cpp**](mpisplitvlasovsolver_8cpp.md)

[Go to the documentation of this file](mpisplitvlasovsolver_8cpp.md)


```C++
// SPDX-License-Identifier: MIT

#include "mpisplitvlasovsolver.hpp"

MpiSplitVlasovSolver::MpiSplitVlasovSolver(
        IAdvectionSpatial<GeometryVxVyXY, GridX> const& advec_x,
        IAdvectionSpatial<GeometryVxVyXY, GridY> const& advec_y,
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

DFieldSpVxVyXY MpiSplitVlasovSolver::operator()(
        DFieldSpVxVyXY const allfdistribu_v2Dsplit,
        DConstFieldXY const electric_field_x,
        DConstFieldXY const electric_field_y,
        double const dt) const
{
    IdxRangeSpVxVyXY idxrange_v2Dsplit(m_transpose.get_local_idx_range<V2DSplit>());
    IdxRangeSpXYVxVy idxrange_x2Dsplit(m_transpose.get_local_idx_range<X2DSplit>());
    IdxRangeXY idx_range_xy_v2Dsplit(idxrange_v2Dsplit);
    DFieldMemSpXYVxVy allfdistribu_x2Dsplit_alloc(idxrange_x2Dsplit);
    DFieldSpXYVxVy allfdistribu_x2Dsplit = get_field(allfdistribu_x2Dsplit_alloc);

    // Create contiguous memory space to contain the relevant section of the electric field
    DFieldMemXY local_electric_field_x(idx_range_xy_v2Dsplit);
    DFieldMemXY local_electric_field_y(idx_range_xy_v2Dsplit);
    ddc::parallel_deepcopy(
            get_field(local_electric_field_x),
            electric_field_x[idx_range_xy_v2Dsplit]);
    ddc::parallel_deepcopy(
            get_field(local_electric_field_y),
            electric_field_y[idx_range_xy_v2Dsplit]);

    // Advect in spatial dimensions
    m_advec_x(allfdistribu_v2Dsplit, dt / 2);
    m_advec_y(allfdistribu_v2Dsplit, dt / 2);
    // Swap to vxvy contiguous layout
    m_transpose(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_x2Dsplit,
            get_const_field(allfdistribu_v2Dsplit));
    // Advect in velocity dimensions
    m_advec_vx(allfdistribu_x2Dsplit, get_const_field(local_electric_field_x), dt / 2);
    m_advec_vy(allfdistribu_x2Dsplit, get_const_field(local_electric_field_y), dt);
    m_advec_vx(allfdistribu_x2Dsplit, get_const_field(local_electric_field_x), dt / 2);
    // Swap to xy contiguous layout
    m_transpose(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_v2Dsplit,
            get_const_field(allfdistribu_x2Dsplit));
    // Advect in spatial dimensions
    m_advec_y(allfdistribu_v2Dsplit, dt / 2);
    m_advec_x(allfdistribu_v2Dsplit, dt / 2);

    return allfdistribu_v2Dsplit;
}
```


