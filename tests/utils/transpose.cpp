// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <ddc_helper.hpp>
#include <transpose.hpp>


namespace {

struct RDimX
{
};
struct RDimY
{
};
struct RDimZ
{
};

using IDimX = ddc::UniformPointSampling<RDimX>;
using IDimY = ddc::UniformPointSampling<RDimY>;
using IDimZ = ddc::UniformPointSampling<RDimZ>;

using IndexX = ddc::DiscreteElement<IDimX>;
using IndexY = ddc::DiscreteElement<IDimY>;
using IndexZ = ddc::DiscreteElement<IDimZ>;
using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;
using IndexYX = ddc::DiscreteElement<IDimY, IDimX>;
using IndexXYZ = ddc::DiscreteElement<IDimX, IDimY, IDimZ>;
using IndexXZY = ddc::DiscreteElement<IDimX, IDimZ, IDimY>;
using IndexZXY = ddc::DiscreteElement<IDimZ, IDimX, IDimY>;
using IndexZYX = ddc::DiscreteElement<IDimZ, IDimY, IDimX>;

using VectX = ddc::DiscreteVector<IDimX>;
using VectY = ddc::DiscreteVector<IDimY>;
using VectZ = ddc::DiscreteVector<IDimZ>;
using VectXY = ddc::DiscreteVector<IDimX, IDimY>;
using VectYX = ddc::DiscreteVector<IDimY, IDimX>;
using VectXYZ = ddc::DiscreteVector<IDimX, IDimY, IDimZ>;
using VectXZY = ddc::DiscreteVector<IDimX, IDimZ, IDimY>;
using VectZXY = ddc::DiscreteVector<IDimZ, IDimX, IDimY>;
using VectZYX = ddc::DiscreteVector<IDimZ, IDimY, IDimX>;

using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainZ = ddc::DiscreteDomain<IDimZ>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;
using IDomainYX = ddc::DiscreteDomain<IDimY, IDimX>;
using IDomainXYZ = ddc::DiscreteDomain<IDimX, IDimY, IDimZ>;
using IDomainXZY = ddc::DiscreteDomain<IDimX, IDimZ, IDimY>;
using IDomainZXY = ddc::DiscreteDomain<IDimZ, IDimX, IDimY>;
using IDomainZYX = ddc::DiscreteDomain<IDimZ, IDimY, IDimX>;

using CoordX = ddc::Coordinate<RDimX>;
using CoordY = ddc::Coordinate<RDimY>;
using CoordZ = ddc::Coordinate<RDimZ>;
using CoordXY = ddc::Coordinate<RDimX, RDimY>;
using CoordXYZ = ddc::Coordinate<RDimX, RDimY, RDimZ>;

using DFieldXY = ddc::Chunk<double, IDomainXY>;
using DFieldYX = ddc::Chunk<double, IDomainYX>;
using DFieldXYZ = ddc::Chunk<double, IDomainXYZ>;
using DFieldXZY = ddc::Chunk<double, IDomainXZY>;
using DFieldZXY = ddc::Chunk<double, IDomainZXY>;
using DFieldZYX = ddc::Chunk<double, IDomainZYX>;

using DSpanXY = ddc::ChunkSpan<double, IDomainXY>;
using DSpanYX = ddc::ChunkSpan<double, IDomainYX>;
using DSpanXYZ = ddc::ChunkSpan<double, IDomainXYZ>;
using DSpanXZY = ddc::ChunkSpan<double, IDomainXZY>;
using DSpanZXY = ddc::ChunkSpan<double, IDomainZXY>;
using DSpanZYX = ddc::ChunkSpan<double, IDomainZYX>;

template <class IDim>
KOKKOS_FUNCTION typename IDim::continuous_element_type get_coordinate(ddc::DiscreteElement<IDim> x)
{
    using RDim = typename IDim::continuous_dimension_type;
    CoordXYZ const origin(-1.0, 10.0, 100.0);
    ddc::DiscreteElement<IDim> const origin_idx(0);
    CoordXYZ step(0.2, 0.1, 2.0);

    return ddc::select<RDim>(origin) + ddc::select<RDim>(step) * (x - origin_idx);
}

TEST(LayoutTransposition, Transpose2D_Host)
{
    IndexXY start_domain_origin(0, 0);
    VectXY start_domain_size(10, 10);
    IndexYX end_domain_origin(0, 0);
    VectYX end_domain_size(10, 10);
    IDomainXY start_domain(start_domain_origin, start_domain_size);
    IDomainYX end_domain(end_domain_origin, end_domain_size);

    DFieldXY start_values_alloc(start_domain);
    DFieldYX end_values_alloc(end_domain);

    DSpanXY start_values = start_values_alloc.span_view();
    DSpanYX end_values = end_values_alloc.span_view();

    ddc::for_each(start_domain, [&](IndexXY ixy) {
        double coord_x = get_coordinate(ddc::select<IDimX>(ixy));
        double coord_y = get_coordinate(ddc::select<IDimY>(ixy));
        start_values(ixy) = coord_x * coord_x - coord_y * coord_y;
    });

    transpose_layout(Kokkos::DefaultHostExecutionSpace(), end_values, start_values.span_cview());

    ddc::for_each(start_domain, [&](IndexXY ixy) {
        IndexX ix(ixy);
        IndexY iy(ixy);
        EXPECT_EQ(start_values(ix, iy), end_values(ix, iy));
    });
}

static void TestTranspose2D_Device()
{
    IndexXY start_domain_origin(0, 0);
    VectXY start_domain_size(10, 10);
    IndexYX end_domain_origin(0, 0);
    VectYX end_domain_size(10, 10);
    IDomainXY start_domain(start_domain_origin, start_domain_size);
    IDomainYX end_domain(end_domain_origin, end_domain_size);

    device_t<DFieldXY> start_values_alloc(start_domain);
    device_t<DFieldYX> end_values_alloc(end_domain);

    device_t<DSpanXY> start_values = start_values_alloc.span_view();
    device_t<DSpanYX> end_values = end_values_alloc.span_view();

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            start_domain,
            KOKKOS_LAMBDA(IndexXY ixy) {
                double coord_x = get_coordinate(ddc::select<IDimX>(ixy));
                double coord_y = get_coordinate(ddc::select<IDimY>(ixy));
                start_values(ixy) = coord_x * coord_x - coord_y * coord_y;
            });

    transpose_layout(Kokkos::DefaultExecutionSpace(), end_values, start_values.span_cview());

    auto start_values_host = ddc::create_mirror_view_and_copy(start_values.span_cview());
    auto end_values_host = ddc::create_mirror_view_and_copy(end_values.span_cview());

    ddc::for_each(start_domain, [&](IndexXY ixy) {
        IndexX ix(ixy);
        IndexY iy(ixy);
        EXPECT_EQ(start_values_host(ix, iy), end_values_host(ix, iy));
    });
}

TEST(LayoutTransposition, Transpose2D_Device)
{
    TestTranspose2D_Device();
}

TEST(LayoutTransposition, BadTranspose2D)
{
    IndexXY start_domain_origin(0, 0);
    VectXY start_domain_size(10, 10);
    IndexYX end_domain_origin(0, 0);
    VectYX end_domain_size(10, 5);
    IDomainXY start_domain(start_domain_origin, start_domain_size);
    IDomainYX end_domain(end_domain_origin, end_domain_size);

    DFieldXY start_values_alloc(start_domain);
    DFieldYX end_values_alloc(end_domain);

    DSpanXY start_values = start_values_alloc.span_view();
    DSpanYX end_values = end_values_alloc.span_view();

    ddc::for_each(start_domain, [&](IndexXY ixy) {
        double coord_x = get_coordinate(ddc::select<IDimX>(ixy));
        double coord_y = get_coordinate(ddc::select<IDimY>(ixy));
        start_values(ixy) = coord_x * coord_x - coord_y * coord_y;
    });

#ifndef NDEBUG
    EXPECT_DEATH(
            transpose_layout(
                    Kokkos::DefaultHostExecutionSpace(),
                    end_values,
                    start_values.span_cview()),
            "Assertion");
#endif
}

TEST(LayoutTransposition, BatchedTranspose2D)
{
    IndexXYZ start_domain_origin(0, 0, 0);
    VectXYZ start_domain_size(10, 10, 10);
    IndexXZY end_domain_origin(0, 0, 0);
    VectXZY end_domain_size(10, 10, 10);
    IDomainXYZ start_domain(start_domain_origin, start_domain_size);
    IDomainXZY end_domain(end_domain_origin, end_domain_size);

    DFieldXYZ start_values_alloc(start_domain);
    DFieldXZY end_values_alloc(end_domain);

    DSpanXYZ start_values = start_values_alloc.span_view();
    DSpanXZY end_values = end_values_alloc.span_view();

    ddc::for_each(start_domain, [&](IndexXYZ ixyz) {
        double coord_x = get_coordinate(ddc::select<IDimX>(ixyz));
        double coord_y = get_coordinate(ddc::select<IDimY>(ixyz));
        double coord_z = get_coordinate(ddc::select<IDimZ>(ixyz));
        start_values(ixyz) = coord_x + coord_y + coord_z;
    });

    transpose_layout(Kokkos::DefaultHostExecutionSpace(), end_values, start_values.span_cview());

    ddc::for_each(start_domain, [&](IndexXYZ ixyz) {
        IndexX ix(ixyz);
        IndexY iy(ixyz);
        IndexZ iz(ixyz);
        EXPECT_EQ(start_values(ix, iy, iz), end_values(ix, iy, iz));
    });
}

TEST(LayoutTransposition, Permutation)
{
    IndexXYZ start_domain_origin(0, 0, 0);
    VectXYZ start_domain_size(10, 10, 10);
    IndexZXY end_domain_origin(0, 0, 0);
    VectZXY end_domain_size(10, 10, 10);
    IDomainXYZ start_domain(start_domain_origin, start_domain_size);
    IDomainZXY end_domain(end_domain_origin, end_domain_size);

    DFieldXYZ start_values_alloc(start_domain);
    DFieldZXY end_values_alloc(end_domain);

    DSpanXYZ start_values = start_values_alloc.span_view();
    DSpanZXY end_values = end_values_alloc.span_view();

    ddc::for_each(start_domain, [&](IndexXYZ ixyz) {
        double coord_x = get_coordinate(ddc::select<IDimX>(ixyz));
        double coord_y = get_coordinate(ddc::select<IDimY>(ixyz));
        double coord_z = get_coordinate(ddc::select<IDimZ>(ixyz));
        start_values(ixyz) = coord_x + coord_y + coord_z;
    });

    transpose_layout(Kokkos::DefaultHostExecutionSpace(), end_values, start_values.span_cview());

    ddc::for_each(start_domain, [&](IndexXYZ ixyz) {
        IndexX ix(ixyz);
        IndexY iy(ixyz);
        IndexZ iz(ixyz);
        EXPECT_EQ(start_values(ix, iy, iz), end_values(ix, iy, iz));
    });
}

TEST(LayoutTransposition, Transpose3D)
{
    IndexXYZ start_domain_origin(0, 0, 0);
    VectXYZ start_domain_size(10, 10, 10);
    IndexZYX end_domain_origin(0, 0, 0);
    VectZYX end_domain_size(10, 10, 10);
    IDomainXYZ start_domain(start_domain_origin, start_domain_size);
    IDomainZYX end_domain(end_domain_origin, end_domain_size);

    DFieldXYZ start_values_alloc(start_domain);
    DFieldZYX end_values_alloc(end_domain);

    DSpanXYZ start_values = start_values_alloc.span_view();
    DSpanZYX end_values = end_values_alloc.span_view();

    ddc::for_each(start_domain, [&](IndexXYZ ixyz) {
        double coord_x = get_coordinate(ddc::select<IDimX>(ixyz));
        double coord_y = get_coordinate(ddc::select<IDimY>(ixyz));
        double coord_z = get_coordinate(ddc::select<IDimZ>(ixyz));
        start_values(ixyz) = coord_x + coord_y + coord_z;
    });

    transpose_layout(Kokkos::DefaultHostExecutionSpace(), end_values, start_values.span_cview());

    ddc::for_each(start_domain, [&](IndexXYZ ixyz) {
        IndexX ix(ixyz);
        IndexY iy(ixyz);
        IndexZ iz(ixyz);
        EXPECT_EQ(start_values(ix, iy, iz), end_values(ix, iy, iz));
    });
}

} // namespace
