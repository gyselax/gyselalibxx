#include <algorithm>
#include <array>
#include <cmath>
#include <iosfwd>
#include <vector>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "cosine_evaluator.hpp"
#include "polynomial_evaluator.hpp"
#include "spline_error_bounds.hpp"

struct DimX
{
    static constexpr bool PERIODIC = false;
};

static constexpr std::size_t s_degree_x = DEGREE_X;

#if defined(BCL_GREVILLE)
static constexpr ddc::BoundCond s_bcl = ddc::BoundCond::GREVILLE;
#elif defined(BCL_HERMITE)
static constexpr ddc::BoundCond s_bcl = ddc::BoundCond::HERMITE;
#endif

#if defined(BCR_GREVILLE)
static constexpr ddc::BoundCond s_bcr = ddc::BoundCond::GREVILLE;
#elif defined(BCR_HERMITE)
static constexpr ddc::BoundCond s_bcr = ddc::BoundCond::HERMITE;
#endif

#if defined(BSPLINES_TYPE_UNIFORM)
using BSplinesX = ddc::UniformBSplines<DimX, s_degree_x>;
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
using BSplinesX = ddc::NonUniformBSplines<DimX, s_degree_x>;
#endif

using GrevillePoints = ddc::GrevilleInterpolationPoints<BSplinesX, s_bcl, s_bcr>;

using IDimX = GrevillePoints::interpolation_mesh_type;

#if defined(EVALUATOR_COSINE)
using evaluator_type = CosineEvaluator::Evaluator<IDimX>;
#elif defined(EVALUATOR_POLYNOMIAL)
using evaluator_type = PolynomialEvaluator::Evaluator<IDimX, s_degree_x>;
#endif

using IndexX = ddc::DiscreteElement<IDimX>;
using DVectX = ddc::DiscreteVector<IDimX>;
using BsplIndexX = ddc::DiscreteElement<BSplinesX>;
using SplineX = ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX>>;
using FieldX = ddc::Chunk<double, ddc::DiscreteDomain<IDimX>>;
using CoordX = ddc::Coordinate<DimX>;

// Checks that when evaluating the spline at interpolation points one
// recovers values that were used to build the spline
TEST(NonPeriodicSplineBuilderTest, Identity)
{
    CoordX constexpr x0(0.);
    CoordX constexpr xN(1.);
    std::size_t constexpr ncells = 100;

    // 1. Create BSplines
    {
#if defined(BSPLINES_TYPE_UNIFORM)
        ddc::init_discrete_space<BSplinesX>(x0, xN, ncells);
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
        DVectX constexpr npoints(ncells + 1);
        std::vector<CoordX> breaks(npoints);
        double dx = (xN - x0) / ncells;
        for (int i(0); i < npoints; ++i) {
            breaks[i] = CoordX(x0 + i * dx);
        }
        ddc::init_discrete_space<BSplinesX>(breaks);
#endif
    }
    ddc::DiscreteDomain<BSplinesX> const dom_bsplines_x(
            ddc::discrete_space<BSplinesX>().full_domain());
    ddc::DiscreteDomain<ddc::Deriv<DimX>> const derivs_domain
            = ddc::DiscreteDomain<ddc::Deriv<DimX>>(
                    ddc::DiscreteElement<ddc::Deriv<DimX>>(1),
                    ddc::DiscreteVector<ddc::Deriv<DimX>>(s_degree_x / 2));

    // 2. Create a Spline represented by a chunk over BSplines
    // The chunk is filled with garbage data, we need to initialize it
    ddc::Chunk coef(dom_bsplines_x, ddc::KokkosAllocator<double, Kokkos::HostSpace>());

    // 3. Create the interpolation domain
    ddc::init_discrete_space<IDimX>(GrevillePoints::get_sampling());
    ddc::DiscreteDomain<IDimX> interpolation_domain(GrevillePoints::get_domain());

    // 4. Create a SplineBuilder over BSplines using some boundary conditions
    ddc::SplineBuilder<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::HostSpace,
            BSplinesX,
            IDimX,
            s_bcl,
            s_bcr,
            ddc::SplineSolver::GINKGO,
            IDimX>
            spline_builder(interpolation_domain);

    // 5. Allocate and fill a chunk over the interpolation domain
    ddc::Chunk yvals(interpolation_domain, ddc::KokkosAllocator<double, Kokkos::HostSpace>());
    evaluator_type evaluator(interpolation_domain);
    evaluator(yvals.span_view());

    int constexpr shift = s_degree_x % 2; // shift = 0 for even order, 1 for odd order
    ddc::Chunk Sderiv_lhs_alloc(derivs_domain, ddc::HostAllocator<double>());
    ddc::ChunkSpan Sderiv_lhs = Sderiv_lhs_alloc.span_view();
    if (s_bcl == ddc::BoundCond::HERMITE) {
        for (ddc::DiscreteElement<ddc::Deriv<DimX>> const ii : derivs_domain) {
            Sderiv_lhs(ii) = evaluator.deriv(x0, ii - derivs_domain.front() + shift);
        }
    }

    ddc::Chunk Sderiv_rhs_alloc(derivs_domain, ddc::HostAllocator<double>());
    ddc::ChunkSpan Sderiv_rhs = Sderiv_rhs_alloc.span_view();
    if (s_bcr == ddc::BoundCond::HERMITE) {
        for (ddc::DiscreteElement<ddc::Deriv<DimX>> const ii : derivs_domain) {
            Sderiv_rhs(ii) = evaluator.deriv(xN, ii - derivs_domain.front() + shift);
        }
    }

// 6. Finally build the spline by filling `coef`
#if defined(BCL_HERMITE)
    auto deriv_l = std::optional(Sderiv_lhs.span_cview());
#else
    decltype(std::optional(Sderiv_lhs.span_cview())) deriv_l = std::nullopt;
#endif

#if defined(BCR_HERMITE)
    auto deriv_r = std::optional(Sderiv_rhs.span_cview());
#else
    decltype(std::optional(Sderiv_rhs.span_cview())) deriv_r = std::nullopt;
#endif

    spline_builder(coef.span_view(), yvals.span_cview(), deriv_l, deriv_r);

    // 7. Create a SplineEvaluator to evaluate the spline at any point in the domain of the BSplines
    ddc::SplineEvaluator<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::HostSpace,
            BSplinesX,
            IDimX,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
            IDimX>
            spline_evaluator(
                    coef.domain(),
                    ddc::NullExtrapolationRule(),
                    ddc::NullExtrapolationRule());

    ddc::Chunk<ddc::Coordinate<DimX>, ddc::DiscreteDomain<IDimX>> coords_eval(interpolation_domain);
    for (IndexX const ix : interpolation_domain) {
        coords_eval(ix) = ddc::coordinate(ix);
    }

    ddc::Chunk spline_eval(interpolation_domain, ddc::KokkosAllocator<double, Kokkos::HostSpace>());
    spline_evaluator(spline_eval.span_view(), coords_eval.span_cview(), coef.span_cview());

    ddc::Chunk spline_eval_deriv(
            interpolation_domain,
            ddc::KokkosAllocator<double, Kokkos::HostSpace>());
    spline_evaluator
            .deriv(spline_eval_deriv.span_view(), coords_eval.span_cview(), coef.span_cview());

    // 8. Checking errors
    double max_norm_error = 0.;
    double max_norm_error_diff = 0.;
    for (IndexX const ix : interpolation_domain) {
        CoordX const x = ddc::coordinate(ix);

        // Compute error
        double const error = spline_eval(ix) - yvals(ix);
        max_norm_error = std::fmax(max_norm_error, std::fabs(error));

        // Compute error
        double const error_deriv = spline_eval_deriv(ix) - evaluator.deriv(x, 1);
        max_norm_error_diff = std::fmax(max_norm_error_diff, std::fabs(error_deriv));
    }
    ddc::Chunk integral(spline_builder.batch_domain(), ddc::HostAllocator<double>());
    spline_evaluator.integrate(integral.span_view(), coef.span_cview());

    double const max_norm_error_integ = std::fabs(
            integral(ddc::DiscreteElement<>()) - evaluator.deriv(xN, -1) + evaluator.deriv(x0, -1));
    double const max_norm = evaluator.max_norm();
    double const max_norm_diff = evaluator.max_norm(1);
    double const max_norm_int = evaluator.max_norm(-1);
    if constexpr (std::is_same_v<
                          evaluator_type,
                          PolynomialEvaluator::Evaluator<IDimX, s_degree_x>>) {
        EXPECT_LE(max_norm_error / max_norm, 1.0e-14);
        EXPECT_LE(max_norm_error_diff / max_norm_diff, 1.0e-12);
        EXPECT_LE(max_norm_error_integ / max_norm_int, 1.0e-14);
    } else {
        SplineErrorBounds<evaluator_type> error_bounds(evaluator);
        const double h = (xN - x0) / ncells;
        EXPECT_LE(
                max_norm_error,
                std::max(error_bounds.error_bound(h, s_degree_x), 1.0e-14 * max_norm));
        EXPECT_LE(
                max_norm_error_diff,
                std::max(error_bounds.error_bound_on_deriv(h, s_degree_x), 1e-12 * max_norm_diff));
        EXPECT_LE(
                max_norm_error_integ,
                std::max(error_bounds.error_bound_on_int(h, s_degree_x), 1.0e-14 * max_norm_int));
    }
}
