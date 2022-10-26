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
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

#include "cosine_evaluator.hpp"
#include "evaluator_2d.hpp"
#include "polynomial_evaluator.hpp"


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

struct DimX
{
    static constexpr bool PERIODIC = false;
};

static constexpr std::size_t s_degree_x = DEGREE_X;

#if UNIFORM == 1
using BSplinesX = UniformBSplines<DimX, s_degree_x>;
#elif UNIFORM == 0
using BSplinesX = NonUniformBSplines<DimX, s_degree_x>;
#endif

using IDimX = SplineBuilder<BSplinesX, s_bcl, s_bcr>::interpolation_mesh_type;
using IndexX = DiscreteElement<IDimX>;
using DVectX = DiscreteVector<IDimX>;
using CoordX = Coordinate<DimX>;

struct DimY
{
    static constexpr bool PERIODIC = false;
};

static constexpr std::size_t s_degree_y = DEGREE_Y;

#if UNIFORM == 1
using BSplinesY = UniformBSplines<DimY, s_degree_y>;
#elif UNIFORM == 0
using BSplinesY = NonUniformBSplines<DimY, s_degree_y>;
#endif

using IDimY = SplineBuilder<BSplinesY, s_bcl, s_bcr>::interpolation_mesh_type;
using IndexY = DiscreteElement<IDimY>;
using DVectY = DiscreteVector<IDimY>;
using CoordY = Coordinate<DimY>;

using IndexXY = DiscreteElement<IDimX, IDimY>;
using BsplIndexXY = DiscreteElement<BSplinesX, BSplinesY>;
using SplineXY = Chunk<double, DiscreteDomain<BSplinesX, BSplinesY>>;
using FieldXY = Chunk<double, DiscreteDomain<IDimX, IDimY>>;
using CoordXY = Coordinate<DimX, DimY>;

using EvaluatorType = Evaluator2D::Evaluator<
        PolynomialEvaluator::Evaluator<IDimX, s_degree_x>,
        PolynomialEvaluator::Evaluator<IDimY, 0>>;

// Checks that when evaluating the spline at interpolation points one
// recovers values that were used to build the spline
TEST(NonPeriodic2DSplineBuilderTest, Identity)
{
    CoordXY constexpr x0(0., 0.);
    CoordXY constexpr xN(1., 1.);
    std::size_t constexpr ncells1 = 100;
    std::size_t constexpr ncells2 = 50;

    // 1. Create BSplines
    {
#if UNIFORM == 1
        init_discrete_space<BSplinesX>(select<DimX>(x0), select<DimX>(xN), ncells1);
#elif UNIFORM == 0
        DVectX constexpr npoints(ncells1 + 1);
        std::vector<CoordX> breaks(npoints);
        double constexpr dx = (get<DimX>(xN) - get<DimX>(x0)) / ncells1;
        for (int i(0); i < npoints; ++i) {
            breaks[i] = CoordX(get<DimX>(x0) + i * dx);
        }
        init_discrete_space<BSplinesX>(breaks);
#endif
    }
    {
#if UNIFORM == 1
        init_discrete_space<BSplinesY>(select<DimY>(x0), select<DimY>(xN), ncells2);
#elif UNIFORM == 0
        DVectY constexpr npoints(ncells2 + 1);
        std::vector<CoordY> breaks(npoints);
        double constexpr dx = (get<DimY>(xN) - get<DimY>(x0)) / ncells2;
        for (int i(0); i < npoints; ++i) {
            breaks[i] = CoordY(get<DimY>(x0) + i * dx);
        }
        init_discrete_space<BSplinesY>(breaks);
#endif
    }

    DiscreteDomain<BSplinesX, BSplinesY> const dom_bsplines_xy(
            BsplIndexXY(0, 0),
            DiscreteVector<BSplinesX, BSplinesY>(
                    discrete_space<BSplinesX>().size(),
                    discrete_space<BSplinesY>().size()));

    // 2. Create a Spline represented by a chunk over BSplines
    // The chunk is filled with garbage data, we need to initialize it
    SplineXY coef(dom_bsplines_xy);

    // 3. Create a SplineBuilder over BSplines using some boundary conditions
    const SplineBuilder2D<BSplinesX, BSplinesY, s_bcl, s_bcr, s_bcl, s_bcr> spline_builder;
    auto interpolation_domain = spline_builder.interpolation_domain();
    auto interpolation_domain_X = spline_builder.interpolation_domain1();
    auto interpolation_domain_Y = spline_builder.interpolation_domain2();

    // 4. Allocate and fill a chunk over the interpolation domain
    FieldXY yvals(interpolation_domain);
    EvaluatorType evaluator;
    evaluator(yvals.span_view());

    int constexpr shift_X = s_degree_x % 2; // shift = 0 for even order, 1 for odd order
    int constexpr shift_Y = s_degree_y % 2; // shift = 0 for even order, 1 for odd order

    int constexpr n_bc_X = s_degree_x / 2;
    int constexpr n_bc_Y = s_degree_y / 2;

    std::vector<double> deriv_xmin_data_X(n_bc_X * interpolation_domain_Y.size());
    std::vector<double> deriv_xmax_data_X(n_bc_X * interpolation_domain_Y.size());
    DSpan2D deriv_xmin(deriv_xmin_data_X.data(), interpolation_domain_Y.size(), n_bc_X);
    DSpan2D deriv_xmax(deriv_xmax_data_X.data(), interpolation_domain_Y.size(), n_bc_X);
    CDSpan2D c_deriv_xmin(deriv_xmin_data_X.data(), interpolation_domain_Y.size(), n_bc_X);
    CDSpan2D c_deriv_xmax(deriv_xmax_data_X.data(), interpolation_domain_Y.size(), n_bc_X);
    auto yiter = interpolation_domain_Y.begin();
    for (std::size_t ii = 0; ii < interpolation_domain_Y.size(); ++ii) {
        for (std::size_t jj = 0; jj < n_bc_X; ++jj) {
            deriv_xmin(ii, jj)
                    = evaluator.deriv(get<DimX>(x0), double(coordinate(*yiter)), jj + shift_X, 0);
            deriv_xmax(ii, jj)
                    = evaluator.deriv(get<DimX>(xN), double(coordinate(*yiter)), jj + shift_X, 0);
        }
        yiter++;
    }

    std::vector<double> deriv_xmin_data_Y(n_bc_Y * interpolation_domain_X.size());
    std::vector<double> deriv_xmax_data_Y(n_bc_Y * interpolation_domain_X.size());
    DSpan2D deriv_ymin(deriv_xmin_data_Y.data(), interpolation_domain_X.size(), n_bc_Y);
    DSpan2D deriv_ymax(deriv_xmax_data_Y.data(), interpolation_domain_X.size(), n_bc_Y);
    CDSpan2D c_deriv_ymin(deriv_xmin_data_Y.data(), interpolation_domain_X.size(), n_bc_Y);
    CDSpan2D c_deriv_ymax(deriv_xmax_data_Y.data(), interpolation_domain_X.size(), n_bc_Y);
    auto xiter = interpolation_domain_X.begin();
    for (std::size_t ii = 0; ii < interpolation_domain_X.size(); ++ii) {
        for (std::size_t jj = 0; jj < n_bc_Y; ++jj) {
            deriv_ymin(ii, jj)
                    = evaluator.deriv(double(coordinate(*xiter)), get<DimY>(x0), 0, jj + shift_Y);
            deriv_ymax(ii, jj)
                    = evaluator.deriv(double(coordinate(*xiter)), get<DimY>(xN), 0, jj + shift_Y);
        }
        xiter++;
    }

    std::vector<double> mixed_derivs_xmin_ymin_data(n_bc_X * n_bc_Y);
    std::vector<double> mixed_derivs_xmin_ymax_data(n_bc_X * n_bc_Y);
    std::vector<double> mixed_derivs_xmax_ymin_data(n_bc_X * n_bc_Y);
    std::vector<double> mixed_derivs_xmax_ymax_data(n_bc_X * n_bc_Y);
    DSpan2D mixed_derivs_xmin_ymin(mixed_derivs_xmin_ymin_data.data(), n_bc_X, n_bc_Y);
    DSpan2D mixed_derivs_xmin_ymax(mixed_derivs_xmin_ymax_data.data(), n_bc_X, n_bc_Y);
    DSpan2D mixed_derivs_xmax_ymin(mixed_derivs_xmax_ymin_data.data(), n_bc_X, n_bc_Y);
    DSpan2D mixed_derivs_xmax_ymax(mixed_derivs_xmax_ymax_data.data(), n_bc_X, n_bc_Y);
    CDSpan2D c_mixed_derivs_xmin_ymin(mixed_derivs_xmin_ymin_data.data(), n_bc_X, n_bc_Y);
    CDSpan2D c_mixed_derivs_xmin_ymax(mixed_derivs_xmin_ymax_data.data(), n_bc_X, n_bc_Y);
    CDSpan2D c_mixed_derivs_xmax_ymin(mixed_derivs_xmax_ymin_data.data(), n_bc_X, n_bc_Y);
    CDSpan2D c_mixed_derivs_xmax_ymax(mixed_derivs_xmax_ymax_data.data(), n_bc_X, n_bc_Y);
    for (std::size_t ii = 0; ii < n_bc_X; ++ii) {
        for (std::size_t jj = 0; jj < n_bc_Y; ++jj) {
            mixed_derivs_xmin_ymin(ii, jj)
                    = evaluator.deriv(get<DimX>(x0), get<DimY>(x0), ii + shift_X, jj + shift_Y);
            mixed_derivs_xmin_ymax(ii, jj)
                    = evaluator.deriv(get<DimX>(x0), get<DimY>(xN), ii + shift_X, jj + shift_Y);
            mixed_derivs_xmax_ymin(ii, jj)
                    = evaluator.deriv(get<DimX>(xN), get<DimY>(x0), ii + shift_X, jj + shift_Y);
            mixed_derivs_xmax_ymax(ii, jj)
                    = evaluator.deriv(get<DimX>(xN), get<DimY>(xN), ii + shift_X, jj + shift_Y);
        }
    }

    const std::optional<CDSpan2D> deriv_l_X(
            s_bcl == BoundCond::HERMITE ? std::optional(c_deriv_xmin) : std::nullopt);
    const std::optional<CDSpan2D> deriv_r_X(
            s_bcr == BoundCond::HERMITE ? std::optional(c_deriv_xmax) : std::nullopt);
    const std::optional<CDSpan2D> deriv_l_Y(
            s_bcl == BoundCond::HERMITE ? std::optional(c_deriv_ymin) : std::nullopt);
    const std::optional<CDSpan2D> deriv_r_Y(
            s_bcr == BoundCond::HERMITE ? std::optional(c_deriv_ymax) : std::nullopt);
    const std::optional<CDSpan2D> md_xmin_ymin(
            s_bcl == BoundCond::HERMITE ? std::optional(c_mixed_derivs_xmin_ymin) : std::nullopt);
    const std::optional<CDSpan2D> md_xmin_ymax(
            s_bcl == BoundCond::HERMITE && s_bcr == BoundCond::HERMITE
                    ? std::optional(c_mixed_derivs_xmin_ymax)
                    : std::nullopt);
    const std::optional<CDSpan2D> md_xmax_ymin(
            s_bcl == BoundCond::HERMITE && s_bcr == BoundCond::HERMITE
                    ? std::optional(c_mixed_derivs_xmax_ymin)
                    : std::nullopt);
    const std::optional<CDSpan2D> md_xmax_ymax(
            s_bcr == BoundCond::HERMITE ? std::optional(c_mixed_derivs_xmax_ymax) : std::nullopt);

    // 5. Finally build the spline by filling `coef`
    spline_builder(
            coef,
            yvals,
            deriv_l_X,
            deriv_r_X,
            deriv_l_Y,
            deriv_r_Y,
            md_xmin_ymin,
            md_xmax_ymin,
            md_xmin_ymax,
            md_xmax_ymax);

    // 6. Create a SplineEvaluator to evaluate the spline at any point in the domain of the BSplines
    const SplineEvaluator2D<BSplinesX, BSplinesY> spline_evaluator(
            g_null_boundary_2d<BSplinesX, BSplinesY>,
            g_null_boundary_2d<BSplinesX, BSplinesY>,
            g_null_boundary_2d<BSplinesX, BSplinesY>,
            g_null_boundary_2d<BSplinesX, BSplinesY>);

    Chunk<CoordXY, DiscreteDomain<IDimX, IDimY>> coords_eval(interpolation_domain);
    for_each(interpolation_domain, [&](IndexXY const ixy) {
        coords_eval(ixy) = CoordXY(coordinate(select<IDimX>(ixy)), coordinate(select<IDimY>(ixy)));
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

    // 7. Checking errors
    double max_norm_error = 0.;
    double max_norm_error_diff1 = 0.;
    double max_norm_error_diff2 = 0.;
    double max_norm_error_diff12 = 0.;
    for_each(interpolation_domain, [&](IndexXY const ixy) {
        IndexX const ix = select<IDimX>(ixy);
        IndexY const iy = select<IDimY>(ixy);
        CoordX const x = coordinate(ix);
        CoordY const y = coordinate(iy);

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
    EXPECT_LE(max_norm_error, 1.0e-12);
    EXPECT_LE(max_norm_error_diff1, 1.0e-3);
    EXPECT_LE(max_norm_error_diff2, 1.0e-3);
    EXPECT_LE(max_norm_error_diff12, 1.0e-3);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::ScopeGuard scope(argc, argv);
    return RUN_ALL_TESTS();
}
