// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "euler.hpp"
#include "rk2.hpp"
#include "rk3.hpp"
#include "mesh_builder.hpp"
#include "r_theta_test_cases.hpp"

namespace {
struct X
{
    static constexpr bool PERIODIC = false;
};
struct Y
{
    static constexpr bool PERIODIC = false;
};

struct R
{
    static constexpr bool PERIODIC = false;
};
struct Theta
{
    static constexpr bool PERIODIC = true;
};

struct GridR : NonUniformGridBase<R>
{
};
struct GridTheta : NonUniformGridBase<Theta>
{
};

static constexpr int BSDegree = 3;
static constexpr ddc::BoundCond SplineRBoundary = ddc::BoundCond::GREVILLE;
static constexpr ddc::BoundCond SplineThetaBoundary = ddc::BoundCond::PERIODIC;


struct BSplinesR : ddc::NonUniformBSplines<R, BSDegree>
{
};
struct BSplinesTheta : ddc::NonUniformBSplines<Theta, BSDegree>
{
};

using SplineRThetaBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultExecutionSpace,
        typename Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        SplineRBoundary, // boundary at r=0
        SplineRBoundary, // boundary at rmax
        SplineThetaBoundary,
        SplineThetaBoundary,
        ddc::SplineSolver::LAPACK>;

using SplineRThetaEvaluator = ddc::SplineEvaluator2D<
        Kokkos::DefaultExecutionSpace,
        typename Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::NullExtrapolationRule, // boundary at r=0
        ddc::ConstantExtrapolationRule<R, Theta>, // boundary at rmax
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>>;

using SplineInterpPointsR
        = ddc::GrevilleInterpolationPoints<BSplinesR, SplineRBoundary, SplineRBoundary>;
using SplineInterpPointsTheta
        = ddc::GrevilleInterpolationPoints<BSplinesTheta, SplineThetaBoundary, SplineThetaBoundary>;

using IdxRangeRTheta = IdxRange<R, Theta>;
}; // namespace

template <
        template <class FieldMem, class DerivFieldMem, class ExecSpace>
        typename TimeStepperBuilder,
        class LogicalToPhysicalMappingType,
        class LogicalToPseudoPhysicalMappingType>
struct PolarAdvectionFixture : public testing::Test
{
    using LogicalToPhysicalMapping = LogicalToPhysicalMappingType;
    using LogicalToPseudoPhysicalMapping = LogicalToPseudoPhysicalMappingType;
    using X_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_x;
    using Y_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_y;
    using TimeStepper = TimeStepperBuilder<
            DField<IdxRangeRTheta>,
            DVectorFieldMem<IdxRangeRTheta, VectorIndexSet<X_adv, Y_adv>>,
            Kokkos::DefaultExecutionSpace>;
};

template <class LogicalToOtherMapping>
LogicalToOtherMapping init_mapping()
{
    using OtherX = typename LogicalToOtherMapping::cartesian_tag_x;
    using OtherY = typename LogicalToOtherMapping::cartesian_tag_y;
    double x0 = 6.2;
    double y0 = 0.8;
    if constexpr (std::is_same_v<
                          LogicalToOtherMapping,
                          CircularToCartesian<R, Theta, OtherX, OtherY>>) {
        return LogicalToOtherMapping(x0, y0);
    } else if constexpr (std::is_same_v<
                                 LogicalToOtherMapping,
                                 CzarnyToCartesian<R, Theta, OtherX, OtherY>>) {
        return LogicalToOtherMapping(0.3, 1.4, y0, y0);
    }
}

using Cases = ::testing::
        Types<std::tuple<Euler, CartesianToCircular<X, Y, R, Theta>, CartesianToCircular<X, Y, R, Theta>>, >;

TYPED_TEST_SUITE(PolarAdvectionFixture, Cases);

TYPED_TEST(PolarAdvectionFixture, Analytical)
{
    using LogicalToPhysicalMapping = typename TestFixture::LogicalToPhysicalMapping;
    using LogicalToPseudoPhysicalMapping = typename TestFixture::LogicalToPseudoPhysicalMapping;
    using TimeStepper = typename TestFixture::TimeStepper;

    Coord<R> const r_min(0.0);
    Coord<R> const r_max(1.0);
    IdxStep<GridR> const nr_cells(15);

    Coord<Theta> const theta_min(0.0);
    Coord<Theta> const theta_max(2.0 * M_PI);
    IdxStep<GridTheta> const ntheta_cells(10);

    ddc::init_discrete_space<BSplinesR>(
            build_random_non_uniform_break_points(r_min, r_max, nr_cells, 0.5));
    ddc::init_discrete_space<BSplinesTheta>(
            build_random_non_uniform_break_points(theta_min, theta_max, ntheta_cells, 0.5));

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::template get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(
            SplineInterpPointsTheta::template get_sampling<GridTheta>());

    IdxRange<GridR> r_idx_range(SplineInterpPointsR::template get_domain<GridR>());
    IdxRange<GridTheta> theta_idx_range(SplineInterpPointsTheta::template get_domain<GridTheta>());
    IdxRange<GridR, GridTheta> idx_range(r_idx_range, theta_idx_range);

    ddc::NullExtrapolationRule r_min_extrap;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrap;
    SplineRThetaBuilder builder(idx_range);
    ddc::ConstantExtrapolationRule<R, Theta> r_max_extrap(r_max);
    SplineRThetaEvaluator evaluator(r_min_extrap, r_max_extrap, theta_extrap, theta_extrap);

    LogicalToPhysicalMapping to_physical = init_mapping<LogicalToPhysicalMapping>();
    LogicalToPseudoPhysicalMapping to_pseudo_physical
            = init_mapping<LogicalToPseudoPhysicalMapping>();

    TimeStepper time_stepper(idx_range);

    SplinePolarFootFinder const
            foot_finder(time_stepper, to_physical, to_pseudo_physical, builder, evaluator);
}
