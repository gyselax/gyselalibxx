// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace {

struct DimX;
using DDimX = ddc::UniformPointSampling<DimX>;

} // namespace

TEST(DiscreteSpace, IsDiscreteSpaceInitialized)
{
    EXPECT_FALSE(ddc::is_discrete_space_initialized<DDimX>());
    ddc::init_discrete_space(DDimX::
                                     init(ddc::Coordinate<DimX>(0),
                                          ddc::Coordinate<DimX>(1),
                                          ddc::DiscreteVector<DDimX>(2)));
    EXPECT_TRUE(ddc::is_discrete_space_initialized<DDimX>());
}
