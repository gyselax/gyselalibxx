#include <cmath>
#include <random>

#include <ddc/ddc.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/polar_bsplines.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

struct R
{
    static constexpr bool PERIODIC = false;
};
struct P
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
static constexpr std::size_t spline_p_degree = DEGREE_P;
static constexpr int continuity = CONTINUITY;

struct BSplinesR : ddc::NonUniformBSplines<R, spline_r_degree>
{
};
#if defined(BSPLINES_TYPE_UNIFORM)
struct BSplinesP : ddc::UniformBSplines<P, spline_p_degree>
{
};
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
struct BSplinesP : ddc::NonUniformBSplines<P, spline_p_degree>
{
};
#endif

using GrevillePointsR = ddc::
        GrevilleInterpolationPoints<BSplinesR, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE>;
using GrevillePointsP = ddc::
        GrevilleInterpolationPoints<BSplinesP, ddc::BoundCond::PERIODIC, ddc::BoundCond::PERIODIC>;

struct GridR : GrevillePointsR::interpolation_discrete_dimension_type
{
};
struct GridP : GrevillePointsP::interpolation_discrete_dimension_type
{
};
struct BSplines : PolarBSplines<BSplinesR, BSplinesP, continuity>
{
};

#if defined(CIRCULAR_MAPPING)
using CircToCart = CircularToCartesian<X, Y, R, P>;
#elif defined(CZARNY_MAPPING)
using CircToCart = CzarnyToCartesian<X, Y, R, P>;
#endif

TEST(PolarSplineTest, ConstantEval)
{
    using PolarCoord = ddc::Coordinate<R, P>;
    using CoordR = ddc::Coordinate<R>;
    using CoordP = ddc::Coordinate<P>;
    using Spline = PolarSpline<BSplines>;
    using Evaluator = PolarSplineEvaluator<BSplines, ddc::NullExtrapolationRule>;
    using BuilderRP = ddc::SplineBuilder2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesP,
            GridR,
            GridP,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            ddc::SplineSolver::LAPACK,
            GridR,
            GridP>;

    using EvaluatorRP = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesP,
            GridR,
            GridP,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
            ddc::PeriodicExtrapolationRule<P>,
            ddc::PeriodicExtrapolationRule<P>,
            GridR,
            GridP>;
    using DiscreteMapping = DiscreteToCartesian<X, Y, BuilderRP, EvaluatorRP>;

    CoordR constexpr r0(0.);
    CoordR constexpr rN(1.);
    CoordP constexpr p0(0.);
    CoordP constexpr pN(2. * M_PI);
    std::size_t constexpr ncells = 20;

    // 1. Create BSplines
    {
        ddc::DiscreteVector<GridR> constexpr npoints_r(ncells + 1);
        std::vector<CoordR> breaks_r(npoints_r);
        const double dr = (rN - r0) / ncells;
        for (int i(0); i < npoints_r; ++i) {
            breaks_r[i] = CoordR(r0 + i * dr);
        }
        ddc::init_discrete_space<BSplinesR>(breaks_r);
#if defined(BSPLINES_TYPE_UNIFORM)
        ddc::init_discrete_space<BSplinesP>(p0, pN, ncells);
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
        ddc::DiscreteVector<GridP> constexpr npoints_p(ncells + 1);
        std::vector<CoordP> breaks_p(npoints_p);
        const double dp = (pN - p0) / ncells;
        for (int i(0); i < npoints_r; ++i) {
            breaks_p[i] = CoordP(p0 + i * dp);
        }
        ddc::init_discrete_space<BSplinesP>(breaks_p);
#endif
    }

    ddc::init_discrete_space<GridR>(GrevillePointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridP>(GrevillePointsP::get_sampling<GridP>());
    ddc::DiscreteDomain<GridR> interpolation_idx_range_R(GrevillePointsR::get_domain<GridR>());
    ddc::DiscreteDomain<GridP> interpolation_idx_range_P(GrevillePointsP::get_domain<GridP>());
    ddc::DiscreteDomain<GridR, GridP>
            interpolation_idx_range(interpolation_idx_range_R, interpolation_idx_range_P);

    BuilderRP builder_rp(interpolation_idx_range);

    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<P> p_extrapolation_rule;
    EvaluatorRP evaluator_rp(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);

#if defined(CIRCULAR_MAPPING)
    CircToCart const coord_changer;
#elif defined(CZARNY_MAPPING)
    CircToCart const coord_changer(0.3, 1.4);
#endif
    DiscreteMapping const mapping
            = DiscreteMapping::analytical_to_discrete(coord_changer, builder_rp, evaluator_rp);
    ddc::init_discrete_space<BSplines>(mapping);

    Spline coef(builder_rp.spline_domain());

    ddc::for_each(coef.singular_spline_coef.domain(), [&](ddc::DiscreteElement<BSplines> const i) {
        coef.singular_spline_coef(i) = 1.0;
    });
    ddc::for_each(
            coef.spline_coef.domain(),
            [&](ddc::DiscreteElement<BSplinesR, BSplinesP> const i) { coef.spline_coef(i) = 1.0; });

    ddc::NullExtrapolationRule extrapolation_rule;
    const Evaluator spline_evaluator(extrapolation_rule);

    std::size_t const n_test_points = 100;
    double const dr = (rN - r0) / n_test_points;
    double const dp = (pN - p0) / n_test_points;

    for (std::size_t i(0); i < n_test_points; ++i) {
        for (std::size_t j(0); j < n_test_points; ++j) {
            PolarCoord const test_point(r0 + i * dr, p0 + j * dp);
            const double val = spline_evaluator(test_point, coef);
            const double deriv_1 = spline_evaluator.deriv_dim_1(test_point, coef);
            const double deriv_2 = spline_evaluator.deriv_dim_2(test_point, coef);

            EXPECT_LE(fabs(val - 1.0), 1.0e-14);
            EXPECT_LE(fabs(deriv_1), 1.0e-13);
            EXPECT_LE(fabs(deriv_2), 1.0e-13);
        }
    }

    Spline integrals(builder_rp.spline_domain());
    ddc::discrete_space<BSplines>().integrals(integrals);
    double area = ddc::transform_reduce(
                          integrals.singular_spline_coef.domain(),
                          0.0,
                          ddc::reducer::sum<double>(),
                          [&](auto idx) { return integrals.singular_spline_coef(idx); })
                  + ddc::transform_reduce(
                          integrals.spline_coef.domain(),
                          0.0,
                          ddc::reducer::sum<double>(),
                          [&](auto idx) { return integrals.spline_coef(idx); });

    EXPECT_NEAR(area, 2 * M_PI, 1e-10);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ::ddc::ScopeGuard ddc_scope(argc, argv);
    return RUN_ALL_TESTS();
}
