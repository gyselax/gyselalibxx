#include <cmath>
#include <random>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "polar_bsplines.hpp"
#include "polar_spline.hpp"
#include "polar_spline_evaluator.hpp"
#include "view.hpp"

namespace {

struct R
{
    static constexpr bool PERIODIC = false;
};
struct Theta
{
    static constexpr bool PERIODIC = true;
};
struct X
{
    static constexpr bool PERIODIC = false;
};
struct Y
{
    static constexpr bool PERIODIC = false;
};

static constexpr std::size_t spline_r_degree = DEGREE_R;
static constexpr std::size_t spline_theta_degree = DEGREE_P;
static constexpr int continuity = CONTINUITY;

struct BSplinesR : ddc::NonUniformBSplines<R, spline_r_degree>
{
};
#if defined(BSPLINES_TYPE_UNIFORM)
struct BSplinesTheta : ddc::UniformBSplines<Theta, spline_theta_degree>
{
};
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
struct BSplinesTheta : ddc::NonUniformBSplines<Theta, spline_theta_degree>
{
};
#endif

using GrevillePointsR = ddc::
        GrevilleInterpolationPoints<BSplinesR, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE>;
using GrevillePointsTheta = ddc::GrevilleInterpolationPoints<
        BSplinesTheta,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;

struct GridR : GrevillePointsR::interpolation_discrete_dimension_type
{
};
struct GridTheta : GrevillePointsTheta::interpolation_discrete_dimension_type
{
};
struct BSplines : PolarBSplines<BSplinesR, BSplinesTheta, continuity>
{
};

#if defined(CIRCULAR_MAPPING)
using CircToCart = CircularToCartesian<R, Theta, X, Y>;
#elif defined(CZARNY_MAPPING)
using CircToCart = CzarnyToCartesian<R, Theta, X, Y>;
#endif

TEST(PolarSplineTest, ConstantEval)
{
    using PolarCoord = Coord<R, Theta>;
    using CoordR = Coord<R>;
    using CoordTheta = Coord<Theta>;
    using SplineMem = PolarSplineMem<BSplines, Kokkos::HostSpace>;
    using Evaluator = PolarSplineEvaluator<BSplines, ddc::NullExtrapolationRule>;
    using BuilderRTheta = ddc::SplineBuilder2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::HostSpace,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            ddc::SplineSolver::LAPACK,
            GridR,
            GridTheta>;

    using EvaluatorRTheta = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::HostSpace,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
            ddc::PeriodicExtrapolationRule<Theta>,
            ddc::PeriodicExtrapolationRule<Theta>,
            GridR,
            GridTheta>;

    CoordR constexpr r0(0.);
    CoordR constexpr rN(1.);
    CoordTheta constexpr theta0(0.);
    CoordTheta constexpr thetaN(2. * M_PI);
    std::size_t constexpr ncells = 20;

    // 1. Create BSplines
    {
        IdxStep<GridR> constexpr npoints_r(ncells + 1);
        std::vector<CoordR> breaks_r(npoints_r);
        const double dr = (rN - r0) / ncells;
        for (int i(0); i < npoints_r; ++i) {
            breaks_r[i] = CoordR(r0 + i * dr);
        }
        ddc::init_discrete_space<BSplinesR>(breaks_r);
#if defined(BSPLINES_TYPE_UNIFORM)
        ddc::init_discrete_space<BSplinesTheta>(theta0, thetaN, ncells);
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
        IdxStep<GridTheta> constexpr npoints_theta(ncells + 1);
        std::vector<CoordTheta> breaks_theta(npoints_theta);
        const double dp = (thetaN - theta0) / ncells;
        for (int i(0); i < npoints_r; ++i) {
            breaks_theta[i] = CoordTheta(theta0 + i * dp);
        }
        ddc::init_discrete_space<BSplinesTheta>(breaks_theta);
#endif
    }

    ddc::init_discrete_space<GridR>(GrevillePointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(GrevillePointsTheta::get_sampling<GridTheta>());
    IdxRange<GridR> interpolation_idx_range_R(GrevillePointsR::get_domain<GridR>());
    IdxRange<GridTheta> interpolation_idx_range_P(GrevillePointsTheta::get_domain<GridTheta>());
    IdxRange<GridR, GridTheta>
            interpolation_idx_range(interpolation_idx_range_R, interpolation_idx_range_P);

    BuilderRTheta builder_rtheta(interpolation_idx_range);

    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    EvaluatorRTheta evaluator_rtheta(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);

#if defined(CIRCULAR_MAPPING)
    CircToCart const coord_changer;
#elif defined(CZARNY_MAPPING)
    CircToCart const coord_changer(0.3, 1.4);
#endif
    DiscreteToCartesianBuilder<X, Y, BuilderRTheta, EvaluatorRTheta> mapping_builder(
            Kokkos::DefaultHostExecutionSpace(),
            coord_changer,
            builder_rtheta,
            evaluator_rtheta);
    DiscreteToCartesian mapping = mapping_builder();
    ddc::init_discrete_space<BSplines>(mapping);

    SplineMem coef(get_spline_idx_range(builder_rtheta));

    ddc::for_each(get_idx_range(coef.singular_spline_coef), [&](Idx<BSplines> const i) {
        coef.singular_spline_coef(i) = 1.0;
    });
    ddc::for_each(get_idx_range(coef.spline_coef), [&](Idx<BSplinesR, BSplinesTheta> const i) {
        coef.spline_coef(i) = 1.0;
    });

    ddc::NullExtrapolationRule extrapolation_rule;
    const Evaluator spline_evaluator(extrapolation_rule);

    std::size_t const n_test_points = 100;
    double const dr = (rN - r0) / n_test_points;
    double const dp = (thetaN - theta0) / n_test_points;

    for (std::size_t i(0); i < n_test_points; ++i) {
        for (std::size_t j(0); j < n_test_points; ++j) {
            PolarCoord const test_point(r0 + i * dr, theta0 + j * dp);
            const double val = spline_evaluator(test_point, get_const_field(coef));
            const double deriv_1 = spline_evaluator.deriv_dim_1(test_point, get_const_field(coef));
            const double deriv_2 = spline_evaluator.deriv_dim_2(test_point, get_const_field(coef));

            EXPECT_LE(fabs(val - 1.0), 1.0e-14);
            EXPECT_LE(fabs(deriv_1), 1.0e-13);
            EXPECT_LE(fabs(deriv_2), 1.0e-13);
        }
    }
}

void test_polar_integrals()
{
    using CoordR = Coord<R>;
    using CoordTheta = Coord<Theta>;
    using SplineMem = PolarSplineMem<BSplines>;
    using Spline = PolarSpline<BSplines>;
    using BuilderRTheta = ddc::SplineBuilder2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            ddc::SplineSolver::LAPACK,
            GridR,
            GridTheta>;

    using EvaluatorRTheta = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
            ddc::PeriodicExtrapolationRule<Theta>,
            ddc::PeriodicExtrapolationRule<Theta>,
            GridR,
            GridTheta>;

    CoordR constexpr r0(0.);
    CoordR constexpr rN(1.);
    CoordTheta constexpr theta0(0.);
    CoordTheta constexpr thetaN(2. * M_PI);
    std::size_t constexpr ncells = 20;

    // 1. Create BSplines
    {
        IdxStep<GridR> constexpr npoints_r(ncells + 1);
        std::vector<CoordR> breaks_r(npoints_r);
        const double dr = (rN - r0) / ncells;
        for (int i(0); i < npoints_r; ++i) {
            breaks_r[i] = CoordR(r0 + i * dr);
        }
        ddc::init_discrete_space<BSplinesR>(breaks_r);
#if defined(BSPLINES_TYPE_UNIFORM)
        ddc::init_discrete_space<BSplinesTheta>(theta0, thetaN, ncells);
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
        IdxStep<GridTheta> constexpr npoints_theta(ncells + 1);
        std::vector<CoordTheta> breaks_theta(npoints_theta);
        const double dp = (thetaN - theta0) / ncells;
        for (int i(0); i < npoints_r; ++i) {
            breaks_theta[i] = CoordTheta(theta0 + i * dp);
        }
        ddc::init_discrete_space<BSplinesTheta>(breaks_theta);
#endif
    }

    ddc::init_discrete_space<GridR>(GrevillePointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(GrevillePointsTheta::get_sampling<GridTheta>());
    IdxRange<GridR> interpolation_idx_range_R(GrevillePointsR::get_domain<GridR>());
    IdxRange<GridTheta> interpolation_idx_range_P(GrevillePointsTheta::get_domain<GridTheta>());
    IdxRange<GridR, GridTheta>
            interpolation_idx_range(interpolation_idx_range_R, interpolation_idx_range_P);

    BuilderRTheta builder_rtheta(interpolation_idx_range);

    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    EvaluatorRTheta evaluator_rtheta(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);

#if defined(CIRCULAR_MAPPING)
    CircToCart const coord_changer;
#elif defined(CZARNY_MAPPING)
    CircToCart const coord_changer(0.3, 1.4);
#endif
    DiscreteToCartesianBuilder<X, Y, BuilderRTheta, EvaluatorRTheta> mapping_builder(
            Kokkos::DefaultHostExecutionSpace(),
            coord_changer,
            builder_rtheta,
            evaluator_rtheta);
    DiscreteToCartesian mapping = mapping_builder();
    ddc::init_discrete_space<BSplines>(mapping);

    SplineMem bspline_integrals_alloc(get_spline_idx_range(builder_rtheta));
    Spline bspline_integrals(bspline_integrals_alloc);
    integrals(Kokkos::DefaultExecutionSpace(), bspline_integrals);
    double area = ddc::parallel_transform_reduce(
                          Kokkos::DefaultExecutionSpace(),
                          get_idx_range(bspline_integrals.singular_spline_coef),
                          0.0,
                          ddc::reducer::sum<double>(),
                          KOKKOS_LAMBDA(Idx<BSplines> idx) {
                              return bspline_integrals.singular_spline_coef(idx);
                          })
                  + ddc::parallel_transform_reduce(
                          Kokkos::DefaultExecutionSpace(),
                          get_idx_range(bspline_integrals.spline_coef),
                          0.0,
                          ddc::reducer::sum<double>(),
                          KOKKOS_LAMBDA(Idx<BSplinesR, BSplinesTheta> idx) {
                              return bspline_integrals.spline_coef(idx);
                          });

    EXPECT_NEAR(area, 2 * M_PI, 1e-10);
}

TEST(PolarSplineTest, Integrals)
{
    test_polar_integrals();
}

} // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ::ddc::ScopeGuard ddc_scope(argc, argv);
    return RUN_ALL_TESTS();
}
