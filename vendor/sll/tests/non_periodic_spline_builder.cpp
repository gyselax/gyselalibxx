#include <algorithm>
#include <array>
#include <cmath>
#include <iosfwd>
#include <vector>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

#include "cosine_evaluator.hpp"
#include "polynomial_evaluator.hpp"

struct DimX
{
    static constexpr bool PERIODIC = false;
};

static constexpr std::size_t s_degree_x = DEGREE_X;

#if BCL == GREVILLE
static constexpr BoundCond s_bcl = BoundCond::GREVILLE;
#elif BCL == HERMITE
static constexpr BoundCond s_bcl = BoundCond::HERMITE;
#endif

#if BCR == GREVILLE
static constexpr BoundCond s_bcr = BoundCond::GREVILLE;
#elif BCR == HERMITE
static constexpr BoundCond s_bcr = BoundCond::HERMITE;
#endif

#if UNIFORM == 1
using BSplinesX = UniformBSplines<DimX, s_degree_x>;
#elif UNIFORM == 0
using BSplinesX = NonUniformBSplines<DimX, s_degree_x>;
#endif

using IDimX = SplineBuilder<BSplinesX, s_bcl, s_bcr>::interpolation_mesh_type;

#if EVALUATOR == Cosine
using evaluator_type = CosineEvaluator::Evaluator<IDimX>;
#elif EVALUATOR == Polynomial
using evaluator_type = PolynomialEvaluator::Evaluator<IDimX>;
#endif

using IndexX = DiscreteElement<IDimX>;
using DVectX = DiscreteVector<IDimX>;
using BsplIndexX = DiscreteElement<BSplinesX>;
using SplineX = Chunk<double, DiscreteDomain<BSplinesX>>;
using FieldX = Chunk<double, DiscreteDomain<IDimX>>;
using CoordX = Coordinate<DimX>;

// Checks that when evaluating the spline at interpolation points one
// recovers values that were used to build the spline
TEST(NonPeriodicSplineBuilderTest, Identity)
{
    CoordX constexpr x0(0.);
    CoordX constexpr xN(1.);
    std::size_t constexpr ncells = 100;

    // 1. Create BSplines
    {
#if UNIFORM == 1
        init_discrete_space<BSplinesX>(x0, xN, ncells);
#elif UNIFORM == 0
        DVectX constexpr npoints(ncells + 1);
        std::vector<CoordX> breaks(npoints);
        double dx = (xN - x0) / ncells;
        for (int i(0); i < npoints; ++i) {
            breaks[i] = CoordX(x0 + i * dx);
        }
        init_discrete_space<BSplinesX>(breaks);
#endif
    }
    DiscreteDomain<BSplinesX> const dom_bsplines_x(discrete_space<BSplinesX>().full_domain());

    // 2. Create a Spline represented by a chunk over BSplines
    // The chunk is filled with garbage data, we need to initialize it
    SplineX coef(dom_bsplines_x);

    // 3. Create a SplineBuilder over BSplines using some boundary conditions
    SplineBuilder<BSplinesX, s_bcl, s_bcr> spline_builder;
    auto interpolation_domain = spline_builder.interpolation_domain();

    // 4. Allocate and fill a chunk over the interpolation domain
    FieldX yvals(interpolation_domain);
    evaluator_type evaluator;
    evaluator(yvals.span_view());

    int constexpr shift = s_degree_x % 2; // shift = 0 for even order, 1 for odd order
    std::array<double, s_degree_x / 2> Sderiv_lhs_data;
    DSpan1D Sderiv_lhs(Sderiv_lhs_data.data(), Sderiv_lhs_data.size());
    std::optional<DSpan1D> deriv_l;
    if (s_bcl == BoundCond::HERMITE) {
        for (std::size_t ii = 0; ii < Sderiv_lhs.extent(0); ++ii) {
            Sderiv_lhs(ii) = evaluator.deriv(x0, ii + shift);
        }
        deriv_l = Sderiv_lhs;
    }

    std::array<double, s_degree_x / 2> Sderiv_rhs_data;
    DSpan1D Sderiv_rhs(Sderiv_rhs_data.data(), Sderiv_rhs_data.size());
    std::optional<DSpan1D> deriv_r;
    if (s_bcr == BoundCond::HERMITE) {
        for (std::size_t ii = 0; ii < Sderiv_rhs.extent(0); ++ii) {
            Sderiv_rhs(ii) = evaluator.deriv(xN, ii + shift);
        }
        deriv_r = Sderiv_rhs;
    }

    // 5. Finally build the spline by filling `coef`
    spline_builder(coef, yvals, deriv_l, deriv_r);

    // 6. Create a SplineEvaluator to evaluate the spline at any point in the domain of the BSplines
    SplineEvaluator<BSplinesX>
            spline_evaluator(g_null_boundary<BSplinesX>, g_null_boundary<BSplinesX>);

    FieldX coords_eval(interpolation_domain);
    for (IndexX const ix : interpolation_domain) {
        coords_eval(ix) = coordinate(ix);
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
        CoordX const x = coordinate(ix);

        // Compute error
        double const error = spline_eval(ix) - yvals(ix);
        max_norm_error = std::fmax(max_norm_error, std::fabs(error));

        // Compute error
        double const error_deriv = spline_eval_deriv(ix) - evaluator.deriv(x, 1);
        max_norm_error_diff = std::fmax(max_norm_error_diff, std::fabs(error_deriv));
    }
    EXPECT_LE(max_norm_error, 1.0e-12);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::ScopeGuard scope(argc, argv);
    return RUN_ALL_TESTS();
}
