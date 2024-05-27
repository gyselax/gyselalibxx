// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "discrete_space.hpp"

struct DimX;
struct DDimX : ddc::UniformPointSampling<DimX>
{
};

void do_not_optimize_away_discrete_space_tests() {}

TEST(DiscreteSpace, Initialization)
{
    ddc::init_discrete_space<DDimX>(DDimX::init<DDimX>(
            ddc::Coordinate<DimX>(0),
            ddc::Coordinate<DimX>(1),
            ddc::DiscreteVector<DDimX>(2)));
}
