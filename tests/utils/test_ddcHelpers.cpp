// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

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
    using Idx = Idx<GridUniform>;
    using IdxStep = IdxStep<GridUniform>;
    using IdxRange = IdxRange<GridUniform>;

    Coord1D x_min(-1.0);
    Coord1D x_max(1.0);
    IdxStep x_size(5);
    double x_len = x_max - x_min;

    ddc::init_discrete_space<GridUniform>(GridUniform::init(x_min, x_max, x_size + 1));

    IdxRange dom(Idx(0), x_size);

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
            Coord1D restricted = ddcHelper::restrict_to_idx_range(test, dom);
            EXPECT_LE(ddc::get<Tag>(restricted), x_max);
            EXPECT_GE(ddc::get<Tag>(restricted), x_min);
            EXPECT_NEAR(ddc::get<Tag>(expected[i]), ddc::get<Tag>(restricted), 1e-8);
        }
    }
}

TEST(DDCHelper, NonUniformPeriodicRestriction)
{
    using Idx = Idx<GridNonUniform>;
    using IdxStep = IdxStep<GridNonUniform>;
    using IdxRange = IdxRange<GridNonUniform>;

    Coord1D x_min(-1.0);
    Coord1D x_max(1.0);
    IdxStep x_size(5);
    double x_len = x_max - x_min;

    std::vector<double> nu_points(6);
    for (int i(0); i < 5; ++i) {
        nu_points[i] = x_min + i * 0.4;
    }
    nu_points[5] = x_max;

    ddc::init_discrete_space<GridNonUniform>(nu_points);

    IdxRange dom(Idx(0), x_size);

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
            Coord1D restricted = ddcHelper::restrict_to_idx_range(test, dom);
            EXPECT_LE(ddc::get<Tag>(restricted), x_max);
            EXPECT_GE(ddc::get<Tag>(restricted), x_min);
            EXPECT_NEAR(ddc::get<Tag>(expected[i]), ddc::get<Tag>(restricted), 1e-8);
        }
    }
}
