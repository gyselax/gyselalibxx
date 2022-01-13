#include <algorithm>
#include <array>
#include <cmath>
#include <iosfwd>
#include <vector>

#include <experimental/mdspan>

#include <ddc/Chunk>
#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/DiscreteDomain>
#include <ddc/NonUniformDiscretization>
#include <ddc/UniformDiscretization>

#include <sll/bsplines_uniform.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

#include "cosine_evaluator.hpp"
#include "test_utils.hpp"

template <class T>
struct PeriodicSplineBuilderTestFixture;

template <std::size_t D, class Evaluator>
struct PeriodicSplineBuilderTestFixture<
        std::tuple<std::integral_constant<std::size_t, D>, Evaluator>> : public testing::Test
{
    struct DimX
    {
        static constexpr bool PERIODIC = true;
    };
    static constexpr std::size_t s_degree = D;
    using IDimX = typename SplineBuilder<
            UniformBSplines<DimX, D>,
            BoundCond::PERIODIC,
            BoundCond::PERIODIC>::interpolation_mesh_type;
    using evaluator_type = typename Evaluator::template Evaluator<IDimX>;
};

using Cases = tuple_to_types_t<cartesian_product_t<
        std::integer_sequence<std::size_t, 1, 2, 3, 4, 5, 6>,
        std::tuple<CosineEvaluator>>>;

TYPED_TEST_SUITE(PeriodicSplineBuilderTestFixture, Cases);

TYPED_TEST(PeriodicSplineBuilderTestFixture, Constructor)
{
    using DimX = typename TestFixture::DimX;
    using BSplinesX = UniformBSplines<DimX, TestFixture::s_degree>;

    init_discretization<BSplinesX>(Coordinate<DimX>(0.), Coordinate<DimX>(0.02), 100);

    SplineBuilder<BSplinesX, BoundCond::PERIODIC, BoundCond::PERIODIC> spline_builder;
    spline_builder.interpolation_domain();
}

// Checks that when evaluating the spline at interpolation points one
// recovers values that were used to build the spline
TYPED_TEST(PeriodicSplineBuilderTestFixture, Identity)
{
    using DimX = typename TestFixture::DimX;
    using IDimX = typename TestFixture::IDimX;
    using IndexX = DiscreteCoordinate<IDimX>;
    using BSplinesX = UniformBSplines<DimX, TestFixture::s_degree>;
    using SplineX = Chunk<double, DiscreteDomain<BSplinesX>>;
    using FieldX = Chunk<double, DiscreteDomain<IDimX>>;
    using CoordX = Coordinate<DimX>;

    CoordX constexpr x0(0.);
    CoordX constexpr xN(1.);
    std::size_t constexpr ncells = 100;
    IndexX constexpr npoints(ncells + 1);

    // 1. Create BSplines
    init_discretization<BSplinesX>(x0, xN, npoints);
    DiscreteDomain<BSplinesX> const dom_bsplines_x(
            DiscreteVector<BSplinesX>(discretization<BSplinesX>().size()));

    // 2. Create a Spline represented by a chunk over BSplines
    // The chunk is filled with garbage data, we need to initialize it
    SplineX coef(dom_bsplines_x);

    // 3. Create a SplineBuilder over BSplines using some boundary conditions
    SplineBuilder<BSplinesX, BoundCond::PERIODIC, BoundCond::PERIODIC> spline_builder;
    auto interpolation_domain = spline_builder.interpolation_domain();

    // 4. Allocate and fill a chunk over the interpolation domain
    FieldX yvals(interpolation_domain);
    typename TestFixture::evaluator_type evaluator;
    evaluator(yvals);

    // 5. Finally build the spline by filling `coef`
    spline_builder(coef, yvals);

    // 6. Create a SplineEvaluator to evaluate the spline at any point in the domain of the BSplines
    SplineEvaluator<BSplinesX> spline_evaluator(NullBoundaryValue::value, NullBoundaryValue::value);

    FieldX coords_eval(interpolation_domain);
    for (IndexX const ix : interpolation_domain) {
        coords_eval(ix) = to_real(ix);
    }

    FieldX spline_eval(interpolation_domain);
    spline_evaluator(spline_eval.span_view(), coords_eval.span_cview(), coef.span_cview());

    FieldX spline_eval_deriv(interpolation_domain);
    spline_evaluator
            .deriv(spline_eval_deriv.span_view(), coords_eval.span_cview(), coef.span_cview());

    // 7. Checking errors
    double max_norm_error = 0.;
    double max_norm_error_diff = 0.;
    for (IndexX const ix : interpolation_domain) {
        CoordX const x = to_real(ix);

        // Compute error
        double const error = spline_eval(ix) - yvals(ix);
        max_norm_error = std::fmax(max_norm_error, std::fabs(error));

        // Compute error
        double const error_deriv = spline_eval_deriv(ix) - evaluator.deriv(x, 1);
        max_norm_error_diff = std::fmax(max_norm_error_diff, std::fabs(error_deriv));
    }
    EXPECT_LE(max_norm_error, 1.0e-12);
}
