#include <cmath>
#include <random>

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/polar_bsplines.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

struct DimR
{
    static constexpr bool PERIODIC = false;
};
struct DimP
{
    static constexpr bool PERIODIC = true;
};
struct DimX
{
    static constexpr bool PERIODIC = false;
};
struct DimY
{
    static constexpr bool PERIODIC = false;
};

static constexpr std::size_t spline_r_degree = DEGREE_R;
static constexpr std::size_t spline_p_degree = DEGREE_P;
static constexpr int continuity = CONTINUITY;

using BSplinesR = NonUniformBSplines<DimR, spline_r_degree>;
#if defined(BSPLINES_TYPE_UNIFORM)
using BSplinesP = UniformBSplines<DimP, spline_p_degree>;
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
using BSplinesP = NonUniformBSplines<DimP, spline_p_degree>;
#endif

using GrevillePointsR
        = GrevilleInterpolationPoints<BSplinesR, BoundCond::GREVILLE, BoundCond::GREVILLE>;
using GrevillePointsP
        = GrevilleInterpolationPoints<BSplinesP, BoundCond::PERIODIC, BoundCond::PERIODIC>;

using IDimR = GrevillePointsR::interpolation_mesh_type;
using IDimP = GrevillePointsP::interpolation_mesh_type;

#if defined(CIRCULAR_MAPPING)
using CircToCart = CircularToCartesian<DimX, DimY, DimR, DimP>;
#elif defined(CZARNY_MAPPING)
using CircToCart = CzarnyToCartesian<DimX, DimY, DimR, DimP>;
#endif

TEST(PolarSplineTest, ConstantEval)
{
    using PolarCoord = ddc::Coordinate<DimR, DimP>;
    using BSplines = PolarBSplines<BSplinesR, BSplinesP, continuity>;
    using CoordR = ddc::Coordinate<DimR>;
    using CoordP = ddc::Coordinate<DimP>;
    using Spline = PolarSpline<BSplines>;
    using Evaluator = PolarSplineEvaluator<BSplines>;
    using BuilderR = SplineBuilder<BSplinesR, IDimR, BoundCond::GREVILLE, BoundCond::GREVILLE>;
    using BuilderP = SplineBuilder<BSplinesP, IDimP, BoundCond::PERIODIC, BoundCond::PERIODIC>;
    using BuilderRP = SplineBuilder2D<BuilderR, BuilderP>;
    using DiscreteMapping = DiscreteToCartesian<DimX, DimY, BuilderRP>;

    CoordR constexpr r0(0.);
    CoordR constexpr rN(1.);
    CoordP constexpr p0(0.);
    CoordP constexpr pN(2. * M_PI);
    std::size_t constexpr ncells = 20;

    // 1. Create BSplines
    {
        ddc::DiscreteVector<IDimR> constexpr npoints_r(ncells + 1);
        std::vector<CoordR> breaks_r(npoints_r);
        const double dr = (rN - r0) / ncells;
        for (int i(0); i < npoints_r; ++i) {
            breaks_r[i] = CoordR(r0 + i * dr);
        }
        ddc::init_discrete_space<BSplinesR>(breaks_r);
#if defined(BSPLINES_TYPE_UNIFORM)
        ddc::init_discrete_space<BSplinesP>(p0, pN, ncells);
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
        ddc::DiscreteVector<IDimP> constexpr npoints_p(ncells + 1);
        std::vector<CoordP> breaks_p(npoints_p);
        const double dp = (pN - p0) / ncells;
        for (int i(0); i < npoints_r; ++i) {
            breaks_p[i] = CoordP(p0 + i * dp);
        }
        ddc::init_discrete_space<BSplinesP>(breaks_p);
#endif
    }

    ddc::init_discrete_space<IDimR>(GrevillePointsR::get_sampling());
    ddc::init_discrete_space<IDimP>(GrevillePointsP::get_sampling());
    ddc::DiscreteDomain<IDimR> interpolation_domain_R(GrevillePointsR::get_domain());
    ddc::DiscreteDomain<IDimP> interpolation_domain_P(GrevillePointsP::get_domain());
    ddc::DiscreteDomain<IDimR, IDimP>
            interpolation_domain(interpolation_domain_R, interpolation_domain_P);

    BuilderR builder_r(interpolation_domain_R);
    BuilderP builder_p(interpolation_domain_P);
    BuilderRP builder_rp(interpolation_domain);

    SplineEvaluator2D<BSplinesR, BSplinesP> evaluator_rp(
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>);

#if defined(CIRCULAR_MAPPING)
    CircToCart const coord_changer;
#elif defined(CZARNY_MAPPING)
    CircToCart const coord_changer(0.3, 1.4);
#endif
    DiscreteMapping const mapping
            = DiscreteMapping::analytical_to_discrete(coord_changer, builder_rp, evaluator_rp);
    ddc::init_discrete_space<BSplines>(mapping, builder_r, builder_p);

    ddc::DiscreteDomain<BSplines> const dom_bsplines_singular(
            ddc::DiscreteElement<BSplines>(0),
            ddc::DiscreteVector<BSplines>(BSplines::n_singular_basis()));
    ddc::DiscreteDomain<BSplinesR> const r_domain(
            ddc::DiscreteElement<BSplinesR>(continuity + 1),
            ddc::DiscreteVector<BSplinesR>(ddc::discrete_space<BSplinesR>().size()));
    ddc::DiscreteDomain<BSplinesP> const p_domain(
            ddc::DiscreteElement<BSplinesP>(0),
            ddc::DiscreteVector<BSplinesP>(ddc::discrete_space<BSplinesP>().size()));
    ddc::DiscreteDomain<BSplinesR, BSplinesP> const dom_2d_bsplines(r_domain, p_domain);

    Spline coef(dom_bsplines_singular, dom_2d_bsplines);

    ddc::for_each(coef.singular_spline_coef.domain(), [&](ddc::DiscreteElement<BSplines> const i) {
        coef.singular_spline_coef(i) = 1.0;
    });
    ddc::for_each(
            coef.spline_coef.domain(),
            [&](ddc::DiscreteElement<BSplinesR, BSplinesP> const i) { coef.spline_coef(i) = 1.0; });

    const Evaluator spline_evaluator(g_polar_null_boundary_2d<BSplines>);

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
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ::ddc::ScopeGuard ddc_scope(argc, argv);
    return RUN_ALL_TESTS();
}
