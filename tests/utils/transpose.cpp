// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "transpose.hpp"


namespace {

struct X
{
};
struct Y
{
};
struct Z
{
};

using GridX = UniformGridBase<X>;
using GridY = UniformGridBase<Y>;
using GridZ = UniformGridBase<Z>;

using IdxX = Idx<GridX>;
using IdxY = Idx<GridY>;
using IdxZ = Idx<GridZ>;
using IdxXY = Idx<GridX, GridY>;
using IdxYX = Idx<GridY, GridX>;
using IdxXYZ = Idx<GridX, GridY, GridZ>;
using IdxXZY = Idx<GridX, GridZ, GridY>;
using IdxZXY = Idx<GridZ, GridX, GridY>;
using IdxZYX = Idx<GridZ, GridY, GridX>;

using IdxStepX = IdxStep<GridX>;
using IdxStepY = IdxStep<GridY>;
using IdxStepZ = IdxStep<GridZ>;
using IdxStepXY = IdxStep<GridX, GridY>;
using IdxStepYX = IdxStep<GridY, GridX>;
using IdxStepXYZ = IdxStep<GridX, GridY, GridZ>;
using IdxStepXZY = IdxStep<GridX, GridZ, GridY>;
using IdxStepZXY = IdxStep<GridZ, GridX, GridY>;
using IdxStepZYX = IdxStep<GridZ, GridY, GridX>;

using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeZ = IdxRange<GridZ>;
using IdxRangeXY = IdxRange<GridX, GridY>;
using IdxRangeYX = IdxRange<GridY, GridX>;
using IdxRangeXYZ = IdxRange<GridX, GridY, GridZ>;
using IdxRangeXZY = IdxRange<GridX, GridZ, GridY>;
using IdxRangeZXY = IdxRange<GridZ, GridX, GridY>;
using IdxRangeZYX = IdxRange<GridZ, GridY, GridX>;

using CoordX = Coord<X>;
using CoordY = Coord<Y>;
using CoordZ = Coord<Z>;
using CoordXY = Coord<X, Y>;
using CoordXYZ = Coord<X, Y, Z>;

using DFieldMemXY = host_t<DFieldMem<IdxRangeXY>>;
using DFieldMemYX = host_t<DFieldMem<IdxRangeYX>>;
using DFieldMemXYZ = host_t<DFieldMem<IdxRangeXYZ>>;
using DFieldMemXZY = host_t<DFieldMem<IdxRangeXZY>>;
using DFieldMemZXY = host_t<DFieldMem<IdxRangeZXY>>;
using DFieldMemZYX = host_t<DFieldMem<IdxRangeZYX>>;

using DFieldXY = host_t<DField<IdxRangeXY>>;
using DFieldYX = host_t<DField<IdxRangeYX>>;
using DFieldXYZ = host_t<DField<IdxRangeXYZ>>;
using DFieldXZY = host_t<DField<IdxRangeXZY>>;
using DFieldZXY = host_t<DField<IdxRangeZXY>>;
using DFieldZYX = host_t<DField<IdxRangeZYX>>;

template <class Grid1D>
KOKKOS_FUNCTION Coord<typename Grid1D::continuous_dimension_type> get_coordinate(Idx<Grid1D> x)
{
    using Dim = typename Grid1D::continuous_dimension_type;
    CoordXYZ const origin(-1.0, 10.0, 100.0);
    Idx<Grid1D> const origin_idx(0);
    CoordXYZ step(0.2, 0.1, 2.0);

    return ddc::select<Dim>(origin) + ddc::select<Dim>(step) * (x - origin_idx);
}

TEST(LayoutTransposition, Transpose2D_Host)
{
    IdxXY start_idx_range_origin(0, 0);
    IdxStepXY start_idx_range_size(10, 10);
    IdxYX end_idx_range_origin(0, 0);
    IdxStepYX end_idx_range_size(10, 10);
    IdxRangeXY start_idx_range(start_idx_range_origin, start_idx_range_size);
    IdxRangeYX end_idx_range(end_idx_range_origin, end_idx_range_size);

    DFieldMemXY start_values_alloc(start_idx_range);
    DFieldMemYX end_values_alloc(end_idx_range);

    DFieldXY start_values = get_field(start_values_alloc);
    DFieldYX end_values = get_field(end_values_alloc);

    ddc::for_each(start_idx_range, [&](IdxXY ixy) {
        double coord_x = get_coordinate(ddc::select<GridX>(ixy));
        double coord_y = get_coordinate(ddc::select<GridY>(ixy));
        start_values(ixy) = coord_x * coord_x - coord_y * coord_y;
    });

    transpose_layout(
            Kokkos::DefaultHostExecutionSpace(),
            end_values,
            get_const_field(start_values));

    ddc::for_each(start_idx_range, [&](IdxXY ixy) {
        IdxX ix(ixy);
        IdxY iy(ixy);
        EXPECT_EQ(start_values(ix, iy), end_values(ix, iy));
    });
}

static void TestTranspose2D_Device()
{
    IdxXY start_idx_range_origin(0, 0);
    IdxStepXY start_idx_range_size(10, 10);
    IdxYX end_idx_range_origin(0, 0);
    IdxStepYX end_idx_range_size(10, 10);
    IdxRangeXY start_idx_range(start_idx_range_origin, start_idx_range_size);
    IdxRangeYX end_idx_range(end_idx_range_origin, end_idx_range_size);

    device_t<DFieldMemXY> start_values_alloc(start_idx_range);
    device_t<DFieldMemYX> end_values_alloc(end_idx_range);

    device_t<DFieldXY> start_values = get_field(start_values_alloc);
    device_t<DFieldYX> end_values = get_field(end_values_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            start_idx_range,
            KOKKOS_LAMBDA(IdxXY ixy) {
                double coord_x = get_coordinate(ddc::select<GridX>(ixy));
                double coord_y = get_coordinate(ddc::select<GridY>(ixy));
                start_values(ixy) = coord_x * coord_x - coord_y * coord_y;
            });

    transpose_layout(Kokkos::DefaultExecutionSpace(), end_values, get_const_field(start_values));

    auto start_values_host = ddc::create_mirror_view_and_copy(get_const_field(start_values));
    auto end_values_host = ddc::create_mirror_view_and_copy(get_const_field(end_values));

    ddc::for_each(start_idx_range, [&](IdxXY ixy) {
        IdxX ix(ixy);
        IdxY iy(ixy);
        EXPECT_EQ(start_values_host(ix, iy), end_values_host(ix, iy));
    });
}

TEST(LayoutTransposition, Transpose2D_Device)
{
    TestTranspose2D_Device();
}

TEST(LayoutTransposition, BadTranspose2D)
{
    IdxXY start_idx_range_origin(0, 0);
    IdxStepXY start_idx_range_size(10, 10);
    IdxYX end_idx_range_origin(0, 0);
    IdxStepYX end_idx_range_size(10, 5);
    IdxRangeXY start_idx_range(start_idx_range_origin, start_idx_range_size);
    IdxRangeYX end_idx_range(end_idx_range_origin, end_idx_range_size);

    DFieldMemXY start_values_alloc(start_idx_range);
    DFieldMemYX end_values_alloc(end_idx_range);

    DFieldXY start_values = get_field(start_values_alloc);
    DFieldYX end_values = get_field(end_values_alloc);

    ddc::for_each(start_idx_range, [&](IdxXY ixy) {
        double coord_x = get_coordinate(ddc::select<GridX>(ixy));
        double coord_y = get_coordinate(ddc::select<GridY>(ixy));
        start_values(ixy) = coord_x * coord_x - coord_y * coord_y;
    });

#ifndef NDEBUG
    EXPECT_DEATH(
            transpose_layout(
                    Kokkos::DefaultHostExecutionSpace(),
                    end_values,
                    get_const_field(start_values)),
            "Assertion");
#endif
}

TEST(LayoutTransposition, BatchedTranspose2D)
{
    IdxXYZ start_idx_range_origin(0, 0, 0);
    IdxStepXYZ start_idx_range_size(10, 10, 10);
    IdxXZY end_idx_range_origin(0, 0, 0);
    IdxStepXZY end_idx_range_size(10, 10, 10);
    IdxRangeXYZ start_idx_range(start_idx_range_origin, start_idx_range_size);
    IdxRangeXZY end_idx_range(end_idx_range_origin, end_idx_range_size);

    DFieldMemXYZ start_values_alloc(start_idx_range);
    DFieldMemXZY end_values_alloc(end_idx_range);

    DFieldXYZ start_values = get_field(start_values_alloc);
    DFieldXZY end_values = get_field(end_values_alloc);

    ddc::for_each(start_idx_range, [&](IdxXYZ ixyz) {
        double coord_x = get_coordinate(ddc::select<GridX>(ixyz));
        double coord_y = get_coordinate(ddc::select<GridY>(ixyz));
        double coord_z = get_coordinate(ddc::select<GridZ>(ixyz));
        start_values(ixyz) = coord_x + coord_y + coord_z;
    });

    transpose_layout(
            Kokkos::DefaultHostExecutionSpace(),
            end_values,
            get_const_field(start_values));

    ddc::for_each(start_idx_range, [&](IdxXYZ ixyz) {
        IdxX ix(ixyz);
        IdxY iy(ixyz);
        IdxZ iz(ixyz);
        EXPECT_EQ(start_values(ix, iy, iz), end_values(ix, iy, iz));
    });
}

TEST(LayoutTransposition, Permutation)
{
    IdxXYZ start_idx_range_origin(0, 0, 0);
    IdxStepXYZ start_idx_range_size(10, 10, 10);
    IdxZXY end_idx_range_origin(0, 0, 0);
    IdxStepZXY end_idx_range_size(10, 10, 10);
    IdxRangeXYZ start_idx_range(start_idx_range_origin, start_idx_range_size);
    IdxRangeZXY end_idx_range(end_idx_range_origin, end_idx_range_size);

    DFieldMemXYZ start_values_alloc(start_idx_range);
    DFieldMemZXY end_values_alloc(end_idx_range);

    DFieldXYZ start_values = get_field(start_values_alloc);
    DFieldZXY end_values = get_field(end_values_alloc);

    ddc::for_each(start_idx_range, [&](IdxXYZ ixyz) {
        double coord_x = get_coordinate(ddc::select<GridX>(ixyz));
        double coord_y = get_coordinate(ddc::select<GridY>(ixyz));
        double coord_z = get_coordinate(ddc::select<GridZ>(ixyz));
        start_values(ixyz) = coord_x + coord_y + coord_z;
    });

    transpose_layout(
            Kokkos::DefaultHostExecutionSpace(),
            end_values,
            get_const_field(start_values));

    ddc::for_each(start_idx_range, [&](IdxXYZ ixyz) {
        IdxX ix(ixyz);
        IdxY iy(ixyz);
        IdxZ iz(ixyz);
        EXPECT_EQ(start_values(ix, iy, iz), end_values(ix, iy, iz));
    });
}

TEST(LayoutTransposition, Transpose3D)
{
    IdxXYZ start_idx_range_origin(0, 0, 0);
    IdxStepXYZ start_idx_range_size(10, 10, 10);
    IdxZYX end_idx_range_origin(0, 0, 0);
    IdxStepZYX end_idx_range_size(10, 10, 10);
    IdxRangeXYZ start_idx_range(start_idx_range_origin, start_idx_range_size);
    IdxRangeZYX end_idx_range(end_idx_range_origin, end_idx_range_size);

    DFieldMemXYZ start_values_alloc(start_idx_range);
    DFieldMemZYX end_values_alloc(end_idx_range);

    DFieldXYZ start_values = get_field(start_values_alloc);
    DFieldZYX end_values = get_field(end_values_alloc);

    ddc::for_each(start_idx_range, [&](IdxXYZ ixyz) {
        double coord_x = get_coordinate(ddc::select<GridX>(ixyz));
        double coord_y = get_coordinate(ddc::select<GridY>(ixyz));
        double coord_z = get_coordinate(ddc::select<GridZ>(ixyz));
        start_values(ixyz) = coord_x + coord_y + coord_z;
    });

    transpose_layout(
            Kokkos::DefaultHostExecutionSpace(),
            end_values,
            get_const_field(start_values));

    ddc::for_each(start_idx_range, [&](IdxXYZ ixyz) {
        IdxX ix(ixyz);
        IdxY iy(ixyz);
        IdxZ iz(ixyz);
        EXPECT_EQ(start_values(ix, iy, iz), end_values(ix, iy, iz));
    });
}

} // namespace
