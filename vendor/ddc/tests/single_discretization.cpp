// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <ddc/experimental/single_discretization.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace ddcexp = ddc::experimental;

namespace {

class DimX;

using CoordX = ddc::Coordinate<DimX>;

struct DDimX : ddcexp::SingleDiscretization<DimX>
{
};

using DElemX = ddc::DiscreteElement<DDimX>;

} // namespace

TEST(SingleDiscretization, ClassSize)
{
    EXPECT_EQ(sizeof(DDimX::Impl<DDimX, Kokkos::HostSpace>), sizeof(double));
}

TEST(SingleDiscretization, Constructor)
{
    CoordX const x(1.);

    ddcexp::SingleDiscretization<DimX>::Impl<DDimX, Kokkos::HostSpace> const ddim_x(x);

    EXPECT_EQ(ddim_x.coordinate(DElemX(0)), x);
}

TEST(SingleDiscretization, Coordinate)
{
    CoordX const x(1.);
    ddc::init_discrete_space<DDimX>(x);
    EXPECT_EQ(ddcexp::coordinate(DElemX(0)), x);
}
