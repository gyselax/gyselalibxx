// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <ddc_helper.hpp>

using namespace ddc;

namespace {

class Tag
{
public:
    static constexpr bool PERIODIC = true;
};

using Coord = Coordinate<Tag>;

struct IDimUniform : UniformPointSampling<Tag>
{
};

struct IDimNonUniform : NonUniformPointSampling<Tag>
{
};

} // namespace

TEST(DDCHelper, UniformPeriodicRestriction)
{
    using Index = DiscreteElement<IDimUniform>;
    using IVect = DiscreteVector<IDimUniform>;
    using IDomain = DiscreteDomain<IDimUniform>;

    Coord x_min(-1.0);
    Coord x_max(1.0);
    IVect x_size(5);
    double x_len = x_max - x_min;

    init_discrete_space<IDimUniform>(IDimUniform::init(x_min, x_max, x_size + 1));

    IDomain dom(Index(0), x_size);

    constexpr int ntest = 10;
    Coord test_start_very_low(x_min - 20 * x_len);
    Coord test_start_below(x_min - 2 * x_len);

    Coord test_start_very_high(x_max + 20 * x_len);
    Coord test_start_above(x_max + 2 * x_len);
    double dx(x_len / ntest);

    std::vector<Coord> expected(ntest);
    for (int i(0); i < ntest; ++i) {
        expected[i] = x_min + dx * i;
    }

    std::vector<Coord> test_cases(
            {test_start_very_low, test_start_below, test_start_above, test_start_very_high});
    for (int j(0); j < 4; ++j) {
        for (int i(0); i < ntest; ++i) {
            Coord test = test_cases[j] + i * dx;
            Coord restricted = ddcHelper::restrict_to_domain(test, dom);
            EXPECT_LE(ddc::get<Tag>(restricted), x_max);
            EXPECT_GE(ddc::get<Tag>(restricted), x_min);
            EXPECT_NEAR(ddc::get<Tag>(expected[i]), ddc::get<Tag>(restricted), 1e-8);
        }
    }
}

TEST(DDCHelper, NonUniformPeriodicRestriction)
{
    using Index = DiscreteElement<IDimNonUniform>;
    using IVect = DiscreteVector<IDimNonUniform>;
    using IDomain = DiscreteDomain<IDimNonUniform>;

    Coord x_min(-1.0);
    Coord x_max(1.0);
    IVect x_size(5);
    double x_len = x_max - x_min;

    std::vector<double> nu_points(6);
    for (int i(0); i < 5; ++i) {
        nu_points[i] = x_min + i * 0.4;
    }
    nu_points[5] = x_max;

    init_discrete_space<IDimNonUniform>(nu_points);

    IDomain dom(Index(0), x_size);

    constexpr int ntest = 10;
    Coord test_start_very_low(x_min - 20 * x_len);
    Coord test_start_below(x_min - 2 * x_len);

    Coord test_start_very_high(x_max + 20 * x_len);
    Coord test_start_above(x_max + 2 * x_len);
    double dx(x_len / ntest);

    std::vector<Coord> expected(ntest);
    for (int i(0); i < ntest; ++i) {
        expected[i] = x_min + dx * i;
    }

    std::vector<Coord> test_cases(
            {test_start_very_low, test_start_below, test_start_above, test_start_very_high});
    for (int j(0); j < 4; ++j) {
        for (int i(0); i < ntest; ++i) {
            Coord test = test_cases[j] + i * dx;
            Coord restricted = ddcHelper::restrict_to_domain(test, dom);
            EXPECT_LE(ddc::get<Tag>(restricted), x_max);
            EXPECT_GE(ddc::get<Tag>(restricted), x_min);
            EXPECT_NEAR(ddc::get<Tag>(expected[i]), ddc::get<Tag>(restricted), 1e-8);
        }
    }
}
