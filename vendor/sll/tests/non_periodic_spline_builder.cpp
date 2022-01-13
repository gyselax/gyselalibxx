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
struct NonPeriodicSplineBuilderTestFixture;

template <std::size_t D, class Evaluator, BoundCond BcL, BoundCond BcR>
struct NonPeriodicSplineBuilderTestFixture<std::tuple<
        std::integral_constant<std::size_t, D>,
        Evaluator,
        std::integral_constant<BoundCond, BcL>,
        std::integral_constant<BoundCond, BcR>>> : public testing::Test
{
    // Needs to be defined here to avoid multiple initializations of the discretization function
    struct DimX
    {
        static constexpr bool PERIODIC = false;
    };
    static constexpr std::size_t s_degree = D;
    static constexpr BoundCond s_bcl = BcL;
    static constexpr BoundCond s_bcr = BcR;
    using IDimX =
            typename SplineBuilder<UniformBSplines<DimX, D>, BcL, BcR>::interpolation_mesh_type;
    using evaluator_type = typename Evaluator::template Evaluator<IDimX>;
};

using Cases = tuple_to_types_t<cartesian_product_t<
        std::integer_sequence<std::size_t, 1, 2, 3, 4, 5, 6>,
        std::tuple<CosineEvaluator>,
        std::integer_sequence<BoundCond, BoundCond::GREVILLE, BoundCond::HERMITE>,
        std::integer_sequence<BoundCond, BoundCond::GREVILLE, BoundCond::HERMITE>>>;

TYPED_TEST_SUITE(NonPeriodicSplineBuilderTestFixture, Cases);

TYPED_TEST(NonPeriodicSplineBuilderTestFixture, Constructor)
{
    using DimX = typename TestFixture::DimX;
    using BSplinesX = UniformBSplines<DimX, TestFixture::s_degree>;

    init_discretization<BSplinesX>(Coordinate<DimX>(0.), Coordinate<DimX>(0.02), 100);

    SplineBuilder<BSplinesX, TestFixture::s_bcl, TestFixture::s_bcr> spline_builder;
    spline_builder.interpolation_domain();
}

// Checks that when evaluating the spline at interpolation points one
// recovers values that were used to build the spline
TYPED_TEST(NonPeriodicSplineBuilderTestFixture, Identity)
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
    SplineBuilder<BSplinesX, TestFixture::s_bcl, TestFixture::s_bcr> spline_builder;
    auto interpolation_domain = spline_builder.interpolation_domain();

    // 4. Allocate and fill a chunk over the interpolation domain
    FieldX yvals(interpolation_domain);
    typename TestFixture::evaluator_type evaluator;
    evaluator(yvals.span_view());

    int constexpr shift = TestFixture::s_degree % 2; // shift = 0 for even order, 1 for odd order
    std::array<double, TestFixture::s_degree / 2> Sderiv_lhs_data;
    DSpan1D Sderiv_lhs(Sderiv_lhs_data.data(), Sderiv_lhs_data.size());
    for (int ii = 0; ii < Sderiv_lhs.extent(0); ++ii) {
        Sderiv_lhs(ii) = evaluator.deriv(x0, ii + shift);
    }
    std::array<double, TestFixture::s_degree / 2> Sderiv_rhs_data;
    DSpan1D Sderiv_rhs(Sderiv_rhs_data.data(), Sderiv_rhs_data.size());
    for (int ii = 0; ii < Sderiv_rhs.extent(0); ++ii) {
        Sderiv_rhs(ii) = evaluator.deriv(xN, ii + shift);
    }
    DSpan1D* deriv_l(TestFixture::s_bcl == BoundCond::HERMITE ? &Sderiv_lhs : nullptr);
    DSpan1D* deriv_r(TestFixture::s_bcr == BoundCond::HERMITE ? &Sderiv_rhs : nullptr);

    // 5. Finally build the spline by filling `coef`
    spline_builder(coef, yvals, deriv_l, deriv_r);

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
