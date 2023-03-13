#include <algorithm>
#include <array>
#include <cmath>
#include <iosfwd>
#include <vector>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_boundary_conditions.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

#include "cosine_evaluator.hpp"
#include "evaluator_2d.hpp"
#include "polynomial_evaluator.hpp"
#include "spline_error_bounds.hpp"

struct DimX
{
    static constexpr bool PERIODIC = true;
};

struct DimY
{
    static constexpr bool PERIODIC = true;
};

static constexpr std::size_t s_degree_x = DEGREE_X;
static constexpr std::size_t s_degree_y = DEGREE_Y;

#if defined(BSPLINES_TYPE_UNIFORM)
using BSplinesX = UniformBSplines<DimX, s_degree_x>;
using BSplinesY = UniformBSplines<DimY, s_degree_y>;
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
using BSplinesX = NonUniformBSplines<DimX, s_degree_x>;
using BSplinesY = NonUniformBSplines<DimY, s_degree_y>;
#endif

using GrevillePointsX
        = GrevilleInterpolationPoints<BSplinesX, BoundCond::PERIODIC, BoundCond::PERIODIC>;
using GrevillePointsY
        = GrevilleInterpolationPoints<BSplinesY, BoundCond::PERIODIC, BoundCond::PERIODIC>;

using IDimX = GrevillePointsX::interpolation_mesh_type;
using IndexX = ddc::DiscreteElement<IDimX>;
using DVectX = ddc::DiscreteVector<IDimX>;
using CoordX = ddc::Coordinate<DimX>;

using IDimY = GrevillePointsY::interpolation_mesh_type;
using IndexY = ddc::DiscreteElement<IDimY>;
using DVectY = ddc::DiscreteVector<IDimY>;
using CoordY = ddc::Coordinate<DimY>;

using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;
using BsplIndexXY = ddc::DiscreteElement<BSplinesX, BSplinesY>;
using SplineXY = ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX, BSplinesY>>;
using FieldXY = ddc::Chunk<double, ddc::DiscreteDomain<IDimX, IDimY>>;
using CoordXY = ddc::Coordinate<DimX, DimY>;

using BuilderX = SplineBuilder<BSplinesX, IDimX, BoundCond::PERIODIC, BoundCond::PERIODIC>;
using BuilderY = SplineBuilder<BSplinesY, IDimY, BoundCond::PERIODIC, BoundCond::PERIODIC>;
using BuilderXY = SplineBuilder2D<BuilderX, BuilderY>;

using EvaluatorType = Evaluator2D::
        Evaluator<CosineEvaluator::Evaluator<IDimX>, CosineEvaluator::Evaluator<IDimY>>;

// Checks that when evaluating the spline at interpolation points one
// recovers values that were used to build the spline
TEST(Periodic2DSplineBuilderTest, Identity)
{
    CoordXY constexpr x0(0., 0.);
    CoordXY constexpr xN(1., 1.);
    std::size_t constexpr ncells1 = 100;
    std::size_t constexpr ncells2 = 50;

    // 1. Create BSplines
    {
#if defined(BSPLINES_TYPE_UNIFORM)
        ddc::init_discrete_space<BSplinesX>(ddc::select<DimX>(x0), ddc::select<DimX>(xN), ncells1);
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
        DVectX constexpr npoints(ncells1 + 1);
        std::vector<CoordX> breaks(npoints);
        double constexpr dx = (ddc::get<DimX>(xN) - ddc::get<DimX>(x0)) / ncells1;
        for (int i(0); i < npoints; ++i) {
            breaks[i] = CoordX(ddc::get<DimX>(x0) + i * dx);
        }
        ddc::init_discrete_space<BSplinesX>(breaks);
#endif
    }
    {
#if defined(BSPLINES_TYPE_UNIFORM)
        ddc::init_discrete_space<BSplinesY>(ddc::select<DimY>(x0), ddc::select<DimY>(xN), ncells2);
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
        DVectY constexpr npoints(ncells2 + 1);
        std::vector<CoordY> breaks(npoints);
        double constexpr dx = (ddc::get<DimY>(xN) - ddc::get<DimY>(x0)) / ncells2;
        for (int i(0); i < npoints; ++i) {
            breaks[i] = CoordY(ddc::get<DimY>(x0) + i * dx);
        }
        ddc::init_discrete_space<BSplinesY>(breaks);
#endif
    }

    ddc::DiscreteDomain<BSplinesX, BSplinesY> const dom_bsplines_xy(
            BsplIndexXY(0, 0),
            ddc::DiscreteVector<BSplinesX, BSplinesY>(
                    ddc::discrete_space<BSplinesX>().size(),
                    ddc::discrete_space<BSplinesY>().size()));

    // 2. Create a Spline represented by a chunk over BSplines
    // The chunk is filled with garbage data, we need to initialize it
    SplineXY coef(dom_bsplines_xy);

    // 3. Create the interpolation domain
    ddc::init_discrete_space<IDimX>(GrevillePointsX::get_sampling());
    ddc::init_discrete_space<IDimY>(GrevillePointsY::get_sampling());
    ddc::DiscreteDomain<IDimX> interpolation_domain_X(GrevillePointsX::get_domain());
    ddc::DiscreteDomain<IDimY> interpolation_domain_Y(GrevillePointsY::get_domain());
    ddc::DiscreteDomain<IDimX, IDimY>
            interpolation_domain(interpolation_domain_X, interpolation_domain_Y);

    // 4. Create a SplineBuilder over BSplines using some boundary conditions
    const BuilderXY spline_builder(interpolation_domain);

    // 5. Allocate and fill a chunk over the interpolation domain
    FieldXY yvals(interpolation_domain);
    EvaluatorType evaluator(interpolation_domain);
    evaluator(yvals.span_view());

    // 6. Finally build the spline by filling `coef`
    spline_builder(coef, yvals);

    // 7. Create a SplineEvaluator to evaluate the spline at any point in the domain of the BSplines
    const SplineEvaluator2D<BSplinesX, BSplinesY> spline_evaluator(
            g_null_boundary_2d<BSplinesX, BSplinesY>,
            g_null_boundary_2d<BSplinesX, BSplinesY>,
            g_null_boundary_2d<BSplinesX, BSplinesY>,
            g_null_boundary_2d<BSplinesX, BSplinesY>);

    ddc::Chunk<CoordXY, ddc::DiscreteDomain<IDimX, IDimY>> coords_eval(interpolation_domain);
    ddc::for_each(interpolation_domain, [&](IndexXY const ixy) {
        coords_eval(ixy) = CoordXY(
                ddc::coordinate(ddc::select<IDimX>(ixy)),
                ddc::coordinate(ddc::select<IDimY>(ixy)));
    });

    FieldXY spline_eval(interpolation_domain);
    spline_evaluator(spline_eval.span_view(), coords_eval.span_cview(), coef.span_cview());

    FieldXY spline_eval_deriv1(interpolation_domain);
    spline_evaluator.deriv_dim_1(
            spline_eval_deriv1.span_view(),
            coords_eval.span_cview(),
            coef.span_cview());

    FieldXY spline_eval_deriv2(interpolation_domain);
    spline_evaluator.deriv_dim_2(
            spline_eval_deriv2.span_view(),
            coords_eval.span_cview(),
            coef.span_cview());

    FieldXY spline_eval_deriv12(interpolation_domain);
    spline_evaluator.deriv_dim_1_and_2(
            spline_eval_deriv12.span_view(),
            coords_eval.span_cview(),
            coef.span_cview());

    // 8. Checking errors
    double max_norm_error = 0.;
    double max_norm_error_diff1 = 0.;
    double max_norm_error_diff2 = 0.;
    double max_norm_error_diff12 = 0.;
    ddc::for_each(interpolation_domain, [&](IndexXY const ixy) {
        IndexX const ix = ddc::select<IDimX>(ixy);
        IndexY const iy = ddc::select<IDimY>(ixy);
        CoordX const x = ddc::coordinate(ix);
        CoordY const y = ddc::coordinate(iy);

        // Compute error
        double const error = spline_eval(ix, iy) - yvals(ix, iy);
        max_norm_error = std::fmax(max_norm_error, std::fabs(error));

        // Compute error
        double const error_deriv1 = spline_eval_deriv1(ix, iy) - evaluator.deriv(x, y, 1, 0);
        max_norm_error_diff1 = std::fmax(max_norm_error_diff1, std::fabs(error_deriv1));

        // Compute error
        double const error_deriv2 = spline_eval_deriv2(ix, iy) - evaluator.deriv(x, y, 0, 1);
        max_norm_error_diff2 = std::fmax(max_norm_error_diff2, std::fabs(error_deriv2));

        // Compute error
        double const error_deriv12 = spline_eval_deriv12(ix, iy) - evaluator.deriv(x, y, 1, 1);
        max_norm_error_diff12 = std::fmax(max_norm_error_diff12, std::fabs(error_deriv12));
    });

    double const max_norm = evaluator.max_norm();
    double const max_norm_diff1 = evaluator.max_norm(1, 0);
    double const max_norm_diff2 = evaluator.max_norm(0, 1);
    double const max_norm_diff12 = evaluator.max_norm(1, 1);

    SplineErrorBounds<EvaluatorType> error_bounds(evaluator);
    const double h1 = (ddc::get<DimX>(xN) - ddc::get<DimX>(x0)) / ncells1;
    const double h2 = (ddc::get<DimY>(xN) - ddc::get<DimY>(x0)) / ncells2;
    EXPECT_LE(
            max_norm_error,
            std::max(error_bounds.error_bound(h1, h2, s_degree_x, s_degree_y), 1.0e-14 * max_norm));
    EXPECT_LE(
            max_norm_error_diff1,
            std::
                    max(error_bounds.error_bound_on_deriv_1(h1, h2, s_degree_x, s_degree_y),
                        1e-12 * max_norm_diff1));
    EXPECT_LE(
            max_norm_error_diff2,
            std::
                    max(error_bounds.error_bound_on_deriv_2(h1, h2, s_degree_x, s_degree_y),
                        1e-12 * max_norm_diff2));
    EXPECT_LE(
            max_norm_error_diff12,
            std::
                    max(error_bounds.error_bound_on_deriv_12(h1, h2, s_degree_x, s_degree_y),
                        1e-10 * max_norm_diff12));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::ddc::ScopeGuard scope(argc, argv);
    return RUN_ALL_TESTS();
}
