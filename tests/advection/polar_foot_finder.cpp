// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>
#include "ddc_aliases.hpp"
#include "r_theta_test_cases.hpp"

class Velocity1DAdvectionTest : public ::testing::Test
{
    struct R;
    struct Theta;

    class GridR : NonUniformGridBase<R>
    {
    };
    class GridTheta : NonUniformGridBase<Theta>
    {
    };
    Velocity1DAdvectionTest()
    {
        Coord<Theta> const theta_min(0.0);
        Coord<Theta> const theta_max(2.0 * M_PI);

        Coord<R> const r_min(0.0);
        Coord<R> const r_max(1.0);
    }
};

TEST(TestSplinePolarFootFinder, Analytical) {}
