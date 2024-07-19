// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <mpilayout.hpp>

namespace {

struct RDimX
{
};
struct RDimY
{
};

struct IDimX : ddc::UniformPointSampling<RDimX>
{
};
struct IDimY : ddc::UniformPointSampling<RDimY>
{
};

using IndexX = ddc::DiscreteElement<IDimX>;
using IndexY = ddc::DiscreteElement<IDimY>;
using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;

using IVectX = ddc::DiscreteVector<IDimX>;
using IVectY = ddc::DiscreteVector<IDimY>;
using IVectXY = ddc::DiscreteVector<IDimX, IDimY>;

using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;
using IDomainYX = ddc::DiscreteDomain<IDimY, IDimX>;

using CoordX = ddc::Coordinate<RDimX>;
using CoordY = ddc::Coordinate<RDimY>;

using XYDistribLayout = MPILayout<IDomainXY, IDimX, IDimY>;
using YDistribLayout = MPILayout<IDomainXY, IDimY>;

TEST(Layout, MinimalDomainDistribution)
{
    IVectX const x_size(8);
    IVectY const y_size(16);

    IndexXY const domain_start(0, 0);
    IVectXY const dom_size(x_size, y_size);
    IDomainXY const global_domain(domain_start, dom_size);

    int n_procs = 4;
    int expected_local_x_extent = x_size.value() / n_procs;

    IndexX const dom_x_start(domain_start);

    XYDistribLayout layout;
    for (int i(0); i < n_procs; ++i) {
        IDomainXY const local_domain = layout.distribute_domain(global_domain, n_procs, i);
        EXPECT_EQ(local_domain.extent<IDimX>().value(), expected_local_x_extent);
        EXPECT_EQ(
                ddc::select<IDimX>(local_domain.front()),
                dom_x_start + i * expected_local_x_extent);
        EXPECT_EQ(local_domain.extent<IDimY>().value(), y_size.value());
    }
}

TEST(Layout, SpreadDomainDistribution)
{
    IVectX const x_size(10);
    IVectY const y_size(15);

    IndexXY const domain_start(0, 0);
    IVectXY const dom_size(x_size, y_size);
    IDomainXY const global_domain(domain_start, dom_size);

    int const n_procs = 6;
    int const expected_procs_x = 2;
    int const expected_procs_y = 3;
    int const expected_local_x_extent = x_size / expected_procs_x;
    int const expected_local_y_extent = y_size / expected_procs_y;

    IndexX const dom_x_start(domain_start);
    IndexY const dom_y_start(domain_start);

    XYDistribLayout layout;
    for (int i(0); i < n_procs; ++i) {
        IDomainXY const local_domain = layout.distribute_domain(global_domain, n_procs, i);
        EXPECT_EQ(local_domain.extent<IDimX>().value(), expected_local_x_extent);
        EXPECT_EQ(local_domain.extent<IDimY>().value(), expected_local_y_extent);
        int const x_rank = i / expected_procs_y;
        int const y_rank = i % expected_procs_y;
        EXPECT_EQ(
                ddc::select<IDimX>(local_domain.front()),
                dom_x_start + x_rank * expected_local_x_extent);
        EXPECT_EQ(
                ddc::select<IDimY>(local_domain.front()),
                dom_y_start + y_rank * expected_local_y_extent);
    }
}

TEST(Layout, DomainSelection)
{
    IVectX const x_size(10);
    IVectY const y_size(15);

    IndexXY const domain_start(0, 0);
    IVectXY const dom_size(x_size, y_size);
    IDomainXY const global_domain(domain_start, dom_size);

    int n_procs = 5;
    int expected_local_y_extent = 3;

    IndexY const dom_y_start(domain_start);

    YDistribLayout layout;
    for (int i(0); i < n_procs; ++i) {
        IDomainXY const local_domain = layout.distribute_domain(global_domain, n_procs, i);
        EXPECT_EQ(local_domain.extent<IDimX>(), global_domain.extent<IDimX>());
        EXPECT_EQ(
                ddc::select<IDimY>(local_domain.front()),
                dom_y_start + i * expected_local_y_extent);
        EXPECT_EQ(local_domain.extent<IDimY>().value(), expected_local_y_extent);
    }
}
} // namespace
