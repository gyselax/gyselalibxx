// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "mesh_builder.hpp"

namespace {

class Tag
{
public:
    static constexpr bool PERIODIC = true;
};

using Coord1D = Coord<Tag>;

struct GridUniform : UniformGridBase<Tag>
{
};

struct GridNonUniform : NonUniformGridBase<Tag>
{
};

} // namespace

TEST(DDCHelper, UniformPeriodicRestriction)
{
    using Grid1D = GridUniform;
    using Idx = Idx<Grid1D>;
    using IdxStep = IdxStep<Grid1D>;
    using IdxRange = IdxRange<Grid1D>;

    Coord1D x_min(-1.0);
    Coord1D x_max(1.0);
    IdxStep npoints(6);
    double x_len = x_max - x_min;

    ddc::init_discrete_space<Grid1D>(Grid1D::init(x_min, x_max, npoints));

    IdxRange idx_range(Idx(0), npoints - 1);

    constexpr int ntest = 10;
    Coord1D test_start_very_low(x_min - 20 * x_len);
    Coord1D test_start_below(x_min - 2 * x_len);

    Coord1D test_start_very_high(x_max + 20 * x_len);
    Coord1D test_start_above(x_max + 2 * x_len);
    double dx(x_len / ntest);

    std::vector<Coord1D> expected(ntest);
    for (int i(0); i < ntest; ++i) {
        expected[i] = x_min + dx * i;
    }

    std::vector<Coord1D> test_cases(
            {test_start_very_low, test_start_below, test_start_above, test_start_very_high});
    for (int j(0); j < 4; ++j) {
        for (int i(0); i < ntest; ++i) {
            Coord1D test = test_cases[j] + i * dx;
            Coord1D restricted = ddcHelper::restrict_to_idx_range(test, idx_range);
            EXPECT_LE(ddc::get<Tag>(restricted), x_max);
            EXPECT_GE(ddc::get<Tag>(restricted), x_min);
            EXPECT_NEAR(ddc::get<Tag>(expected[i]), ddc::get<Tag>(restricted), 1e-8);
        }
    }
}

TEST(DDCHelper, NonUniformPeriodicRestriction)
{
    using Grid1D = GridNonUniform;
    using Idx = Idx<Grid1D>;
    using IdxStep = IdxStep<Grid1D>;
    using IdxRange = IdxRange<Grid1D>;

    Coord1D x_min(-1.0);
    Coord1D x_max(1.0);
    IdxStep ncells(5);
    IdxStep npoints(ncells + 1);
    double x_len = x_max - x_min;

    ddc::init_discrete_space<Grid1D>(build_random_non_uniform_break_points(x_min, x_max, ncells));

    IdxRange idx_range(Idx(0), npoints - 1);

    constexpr int ntest = 10;
    Coord1D test_start_very_low(x_min - 20 * x_len);
    Coord1D test_start_below(x_min - 2 * x_len);

    Coord1D test_start_very_high(x_max + 20 * x_len);
    Coord1D test_start_above(x_max + 2 * x_len);
    double dx(x_len / ntest);

    std::vector<Coord1D> expected(ntest);
    for (int i(0); i < ntest; ++i) {
        expected[i] = x_min + dx * i;
    }

    std::vector<Coord1D> test_cases(
            {test_start_very_low, test_start_below, test_start_above, test_start_very_high});
    for (int j(0); j < 4; ++j) {
        for (int i(0); i < ntest; ++i) {
            Coord1D test = test_cases[j] + i * dx;
            Coord1D restricted = ddcHelper::restrict_to_idx_range(test, idx_range);
            EXPECT_LE(ddc::get<Tag>(restricted), x_max);
            EXPECT_GE(ddc::get<Tag>(restricted), x_min);
            EXPECT_NEAR(ddc::get<Tag>(expected[i]), ddc::get<Tag>(restricted), 1e-8);
        }
    }
}

/**
 * A test for the maximum_distance_between_adjacent_points function,
 * with a uniform grid.
 */
TEST(DDCHelper, ComputeMaxDistanceUniformGrid)
{
    using Grid1D = GridUniform;
    using Idx = Idx<Grid1D>;
    using IdxStep = IdxStep<Grid1D>;
    using IdxRange = IdxRange<Grid1D>;

    Coord1D x_min(0.);
    Coord1D x_max(1.0);
    IdxStep ncells(10);
    IdxStep npoints(ncells + 1);
    double x_len = x_max - x_min;

    ddc::init_discrete_space<Grid1D>(Grid1D::init(x_min, x_max, npoints));

    IdxRange idx_range(Idx(0), npoints);
    double const max_dx(x_len / (ncells));
    double const max_dx_ddchelper(
            ddcHelper::maximum_distance_between_adjacent_points<Grid1D>(idx_range));

    EXPECT_NEAR(max_dx, max_dx_ddchelper, 1e-12);
}

/**
 * A test for the maximum_distance_between_adjacent_points function,
 * with a non uniform grid.
 */
TEST(DDCHelper, ComputeMaxDistanceNonUniformGridFirst)
{
    using Grid1D = GridNonUniform;
    using Idx = Idx<Grid1D>;
    using IdxStep = IdxStep<Grid1D>;
    using IdxRange = IdxRange<Grid1D>;


    std::vector<double> points_list = {0., 0.9, 1., 1.1, 1.5, 2.};
    ddc::init_discrete_space<Grid1D>(points_list);

    IdxStep ncells(5);
    IdxStep npoints(ncells + 1);
    IdxRange idx_range(Idx(0), npoints);

    double const max_dx(0.9);
    double const max_dx_ddchelper(
            ddcHelper::maximum_distance_between_adjacent_points<Grid1D>(idx_range));

    EXPECT_NEAR(max_dx, max_dx_ddchelper, 1e-12);
}

/**
 * A test for the maximum_distance_between_adjacent_points function,
 * with a non uniform grid.
 */
TEST(DDCHelper, ComputeMaxDistanceNonUniformGridMiddle)
{
    using Grid1D = GridNonUniform;
    using Idx = Idx<Grid1D>;
    using IdxStep = IdxStep<Grid1D>;
    using IdxRange = IdxRange<Grid1D>;


    std::vector<double> points_list = {0., 0.1, 1., 1.1, 1.5, 2.};
    ddc::init_discrete_space<Grid1D>(points_list);

    IdxStep ncells(5);
    IdxStep npoints(ncells + 1);
    IdxRange idx_range(Idx(0), npoints);

    double const max_dx(0.9);
    double const max_dx_ddchelper(
            ddcHelper::maximum_distance_between_adjacent_points<Grid1D>(idx_range));

    EXPECT_NEAR(max_dx, max_dx_ddchelper, 1e-12);
}

/**
 * A test for the maximum_distance_between_adjacent_points function,
 * with a non uniform grid.
 */
TEST(DDCHelper, ComputeMaxDistanceNonUniformGridLast)
{
    using Grid1D = GridNonUniform;
    using Idx = Idx<Grid1D>;
    using IdxStep = IdxStep<Grid1D>;
    using IdxRange = IdxRange<Grid1D>;


    std::vector<double> points_list = {0., 0.3, 0.9, 1., 1.1, 2.};
    ddc::init_discrete_space<Grid1D>(points_list);

    int ncells(5);
    IdxStep npoints(ncells + 1);
    IdxRange idx_range(Idx(0), npoints);

    double const max_dx(0.9);
    double const max_dx_ddchelper(
            ddcHelper::maximum_distance_between_adjacent_points<Grid1D>(idx_range));

    EXPECT_NEAR(max_dx, max_dx_ddchelper, 1e-12);
}
