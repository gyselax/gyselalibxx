// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <array>
#include <cmath>
#include <iosfwd>
#include <vector>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

struct DimX
{
    static constexpr bool PERIODIC = true;
};

static constexpr std::size_t s_degree_x = DEGREE_X;

struct BSplinesX : ddc::NonUniformBSplines<DimX, s_degree_x>
{
};

using GrevillePoints = ddc::
        GrevilleInterpolationPoints<BSplinesX, ddc::BoundCond::PERIODIC, ddc::BoundCond::PERIODIC>;

struct IDimX : GrevillePoints::interpolation_mesh_type
{
};

using IndexX = ddc::DiscreteElement<IDimX>;
using DVectX = ddc::DiscreteVector<IDimX>;
using BsplIndexX = ddc::DiscreteElement<BSplinesX>;
using SplineX = ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX>>;
using FieldX = ddc::Chunk<double, ddc::DiscreteDomain<IDimX>>;
using CoordX = ddc::Coordinate<DimX>;

TEST(PeriodicSplineBuilderOrderTest, OrderedPoints)
{
    std::size_t constexpr ncells = 10;

    // 1. Create BSplines
    int constexpr npoints(ncells + 1);
    std::vector<double> d_breaks({0, 0.01, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0});
    std::vector<CoordX> breaks(npoints);
    for (std::size_t i(0); i < npoints; ++i) {
        breaks[i] = CoordX(d_breaks[i]);
    }
    ddc::init_discrete_space<BSplinesX>(breaks);

    // 2. Create the interpolation domain
    ddc::init_discrete_space<IDimX>(GrevillePoints::get_sampling<IDimX>());
    ddc::DiscreteDomain<IDimX> interpolation_domain(GrevillePoints::get_domain<IDimX>());

    double last(ddc::coordinate(interpolation_domain.front()));
    double current;
    for (IndexX const ix : interpolation_domain) {
        current = ddc::coordinate(ix);
        EXPECT_LE(current, ddc::discrete_space<BSplinesX>().rmax());
        EXPECT_GE(current, ddc::discrete_space<BSplinesX>().rmin());
        EXPECT_LE(last, current);
        last = current;
    }
}
