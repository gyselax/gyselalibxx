// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>

#include "ddc_aliases.hpp"
#include "mesh_builder.hpp"
#include "r_theta_test_cases.hpp"

namespace {
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
    using SplineInterpPointsTheta = ddc::
            GrevilleInterpolationPoints<BSplinesTheta, SplineThetaBoundary, SplineThetaBoundary>;
};

template<class TimeStepperType>
class PolarAdvectionFixture : public testing::Test
{
public:
    using TimeStepper = TimeStepperType;

private:
    ddc::NullExtrapolationRule m_r_min_extrap;
    ddc::PeriodicExtrapolationRule<Theta> m_theta_extrap;
    std::unique_ptr<ddc::ConstantExtrapolationRule<R, Theta>> m_r_max_extrap;
    std::unique_ptr<SplineRThetaBuilder> m_builder;
    std::unique_ptr<SplineRThetaEvaluator> m_evaluator;

public:
    PolarAdvectionTest()
    {
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
        IdxRange<GridTheta> theta_idx_range(
                SplineInterpPointsTheta::template get_domain<GridTheta>());
        IdxRange<GridR, GridTheta> idx_range(r_idx_range, theta_idx_range);

        m_builder = std::make_unique<SplineRThetaBuilder>(idx_range);
        m_r_max_extrap = std::make_unique<ddc::ConstantExtrapolationRule<R, Theta>>(r_max);
        m_evaluator = std::make_unique<SplineRThetaEvaluator>(
                m_r_min_extrap,
                *m_r_max_extrap,
                m_theta_extrap,
                m_theta_extrap);
    }
};

using Cases = ::testing::Types<RK3<FieldMem,
        VectorFieldMem,
        Kokkos::DefaultExecutionSpace>>>;

TYPED_TEST_SUITE(PolarBsplineFixture, Cases);

TEST_P(TestSplinePolarFootFinder, Analytical) {}

INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        PolarAdvectionTest,
        testing::Combine(testing::Values<std::size_t>(64), testing::Values<std::size_t>(128)));

