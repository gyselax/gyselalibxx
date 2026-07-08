// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "test_utils.hpp"

template <class T>
struct NonUniformInterpolationPointsFixture;

template <ddc::SplineBuilderClosure BcMin, ddc::SplineBuilderClosure BcMax>
struct NonUniformInterpolationPointsFixture<std::tuple<
        std::integral_constant<ddc::SplineBuilderClosure, BcMin>,
        std::integral_constant<ddc::SplineBuilderClosure, BcMax>>> : public testing::Test
{
    struct X
    {
        static bool constexpr PERIODIC = (BcMin == ddc::SplineBuilderClosure::PERIODIC);
    };
    struct GridX : NonUniformGridBase<X>
    {
    };
    struct BSplinesX : ddc::NonUniformBSplines<X, 3>
    {
    };

    static constexpr ddc::SplineBuilderClosure bc_min = BcMin;
    static constexpr ddc::SplineBuilderClosure bc_max = BcMax;
};

using Cases = tuple_to_types_t<std::tuple<
        std::tuple<
                std::integral_constant<
                        ddc::SplineBuilderClosure,
                        ddc::SplineBuilderClosure::HERMITE>,
                std::integral_constant<
                        ddc::SplineBuilderClosure,
                        ddc::SplineBuilderClosure::HERMITE>>,
        std::tuple<
                std::integral_constant<
                        ddc::SplineBuilderClosure,
                        ddc::SplineBuilderClosure::GREVILLE>,
                std::integral_constant<
                        ddc::SplineBuilderClosure,
                        ddc::SplineBuilderClosure::GREVILLE>>,
        std::tuple<
                std::integral_constant<
                        ddc::SplineBuilderClosure,
                        ddc::SplineBuilderClosure::HERMITE>,
                std::integral_constant<
                        ddc::SplineBuilderClosure,
                        ddc::SplineBuilderClosure::GREVILLE>>,
        std::tuple<
                std::integral_constant<
                        ddc::SplineBuilderClosure,
                        ddc::SplineBuilderClosure::PERIODIC>,
                std::integral_constant<
                        ddc::SplineBuilderClosure,
                        ddc::SplineBuilderClosure::PERIODIC>>>>;

TYPED_TEST_SUITE(NonUniformInterpolationPointsFixture, Cases);

TYPED_TEST(NonUniformInterpolationPointsFixture, CoordinatesMatchAfterInit)
{
    using X = typename TestFixture::X;
    using GridX = typename TestFixture::GridX;
    using BSplinesX = typename TestFixture::BSplinesX;
    using SplineInterpPoints = ddcHelper::
            NonUniformInterpolationPoints<BSplinesX, TestFixture::bc_min, TestFixture::bc_max>;

    const Coord<X> x_min(0.0);
    const Coord<X> x_max(1.0);
    const int n_cells = 4;

    std::vector<Coord<X>> break_points(n_cells + 1);
    for (int i = 0; i <= n_cells; ++i) {
        break_points[i] = x_min + i * (x_max - x_min) / n_cells;
    }

    ddc::init_discrete_space<BSplinesX>(break_points);

    int const npoints = ddc::discrete_space<BSplinesX>().nbasis() - SplineInterpPoints::N_BE_MIN
                        - SplineInterpPoints::N_BE_MAX;

    // Generate npoints evenly distributed points interior to the domain,
    // with at most one per sub-interval, to satisfy the cell constraint.
    std::vector<Coord<X>> interp_points(npoints);
    for (int i = 0; i < npoints; ++i) {
        interp_points[i] = x_min + (i + 0.5) * (x_max - x_min) / npoints;
    }

    ddc::init_discrete_space<GridX>(
            SplineInterpPoints::template get_sampling<GridX>(interp_points));

    IdxRange<GridX> idx_range = SplineInterpPoints::template get_domain<GridX>();
    EXPECT_EQ(idx_range.size(), npoints);

    for (int i(0); i < idx_range.size(); ++i) {
        EXPECT_DOUBLE_EQ(double(ddc::coordinate(idx_range.front() + i)), double(interp_points[i]));
    }

    if constexpr (TestFixture::bc_min == ddc::SplineBuilderClosure::PERIODIC) {
        EXPECT_DOUBLE_EQ(ddcHelper::total_interval_length(idx_range), x_max - x_min);
    }
}
