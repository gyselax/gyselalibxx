// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines/deriv.hpp>

#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "view.hpp"

namespace {
struct X
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = X;
};

struct Y
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Y;
};

struct GridX : UniformGridBase<X>
{
};

struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};

struct GridY : UniformGridBase<Y>
{
};

struct BSplinesY : ddc::UniformBSplines<Y, 3>
{
};

using dX = ddc::Deriv<GridX>;
using DerivX = ddc::Deriv<X>;
using IdxX = Idx<GridX>;
using IdxStepX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;

using dY = ddc::Deriv<GridY>;
using DerivY = ddc::Deriv<Y>;
using IdxY = Idx<GridY>;
using IdxStepY = IdxStep<GridY>;
using IdxRangeY = IdxRange<GridY>;


using IdxXY = Idx<GridX, GridY>;
using IdxStepXY = IdxStep<GridX, GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;

using IdxdXX = Idx<dX, GridX>;
using IdxdYX = Idx<dY, GridX>;
using IdxdYY = Idx<dY, GridY>;
using IdxdXY = Idx<dX, GridY>;

using IdxdXdYXY = Idx<dX, dY, GridX, GridY>;
using IdxRangedXdYXY = IdxRange<dX, dY, GridX, GridY>;

using SplineInterpPointsX = ddc::
        KnotsAsInterpolationPoints<BSplinesX, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;
using SplineInterpPointsY = ddc::
        KnotsAsInterpolationPoints<BSplinesY, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;

using DerivFieldMemdXdYXY = DerivFieldMem<double, IdxRangedXdYXY, 2>;
using DerivFielddXdYXY = DerivField<double, IdxRangedXdYXY>;

} // namespace



TEST(DerivFieldTest, SplineBuilderUse)
{
    // Define meshes -----------------------------------------------------------------------------
    const Coord<X> xmin(0);
    const Coord<X> xmax(1);
    const IdxStepX x_ncells(10);

    const Coord<Y> ymin(0);
    const Coord<Y> ymax(1);
    const IdxStepY y_ncells(10);

    ddc::init_discrete_space<BSplinesX>(xmin, xmax, x_ncells);
    ddc::init_discrete_space<BSplinesY>(ymin, ymax, y_ncells);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridY>(SplineInterpPointsY::get_sampling<GridY>());

    // For the number of derivatives (first, second ... derivatives)
    Idx<dX> idx0_dx(0);
    IdxStep<dX> nelems_dx(1);
    IdxRange<dX> idx_range_dx(idx0_dx, nelems_dx);

    Idx<dY> idx0_dy(0);
    IdxStep<dY> nelems_dy(1);
    IdxRange<dY> idx_range_dy(idx0_dy, nelems_dy);

    // For the point grids
    IdxRangeX idx_range_x(SplineInterpPointsX::template get_domain<GridX>());
    IdxRangeY idx_range_y(SplineInterpPointsY::template get_domain<GridY>());
    IdxRangeXY idx_range_xy(idx_range_x, idx_range_y);
    IdxRangedXdYXY idx_range_dxdy_xy(idx_range_dx, idx_range_dy, idx_range_x, idx_range_y);

    // A subset of the x-derivatives to be retrieved with get_mdspan
    IdxRange<dX> deriv_block_x(Idx<dX>(1), IdxStep<dX>(2));
    // A subset of the y-derivatives to be retrieved with get_mdspan
    IdxRange<dY> deriv_block_y(Idx<dY>(1), IdxStep<dY>(2));
    // A subset of the x and y derivatives to be retrieved with get_mdspan
    IdxRange<dX, dY> deriv_block_xy(deriv_block_x, deriv_block_y);


    // Index range slices where derivatives are defined on the point grids
    IdxRangeSlice<GridX>
            idx_range_slice_dx(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());
    IdxRangeSlice<GridY>
            idx_range_slice_dy(idx_range_y.front(), IdxStepY(2), idx_range_y.extents());

    // Define a DerivField -----------------------------------------------------------------------
    DerivFieldMemdXdYXY
            function_and_derivs_alloc(idx_range_xy, idx_range_slice_dx, idx_range_slice_dy);
    DerivFielddXdYXY function_and_derivs(function_and_derivs_alloc);

    // Extract the function values and its derivatives -------------------------------------------
    Idx<dX> first_dx(1);
    Idx<dY> first_dy(1);

    DField<IdxRange<GridX, GridY>, Kokkos::HostSpace, Kokkos::layout_stride> function_extracted
            = function_and_derivs.get_values_field();

    DField<IdxRange<GridY>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_xmin_extracted
            = function_and_derivs[IdxdXX(first_dx, idx_range_slice_dx.front())];
    DField<IdxRange<GridY>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_xmax_extracted
            = function_and_derivs[IdxdXX(first_dx, idx_range_slice_dx.back())];

    DField<IdxRange<GridX>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_ymin_extracted
            = function_and_derivs[IdxdYY(first_dy, idx_range_slice_dy.front())];
    DField<IdxRange<GridX>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_ymax_extracted
            = function_and_derivs[IdxdYY(first_dy, idx_range_slice_dy.back())];

    /*
        EXAMPLE OF USE FOR THE FIRST DERIVATIVES. 

    // Defined on (dX, GridX, slice GridY)
    detail::ViewNDMaker<3, double, false>::type derivs_x_extracted
            = function_and_derivs.get_mdspan(deriv_block_x);
    // Defined on (dY, slice GridX, GridY)
    detail::ViewNDMaker<3, double, false>::type derivs_y_extracted
            = function_and_derivs.get_mdspan(deriv_block_y);
    */

    // Defined on (dX, dY, slice GridX, GridY)
    detail::ViewNDMaker<4, double, false>::type derivs_xy_extracted
            = function_and_derivs.get_mdspan(deriv_block_xy);


    // Initialise the function and its derivatives -----------------------------------------------
    ddc::for_each(idx_range_xy, [&](IdxXY idx_xy) {
        double const x = ddc::select<X>(ddc::coordinate(idx_xy));
        double const y = ddc::select<Y>(ddc::coordinate(idx_xy));
        function_extracted(idx_xy) = x * x + y * y + 2 * x * y;
    });

    ddc::for_each(idx_range_y, [&](IdxY idx_y) {
        double const y = ddc::coordinate(idx_y);
        double const x_min(xmin);
        double const x_max(xmax);
        derivs_xmin_extracted(idx_y) = 2 * x_min + 2 * y;
        derivs_xmax_extracted(idx_y) = 2 * x_max + 2 * y;
    });

    ddc::for_each(idx_range_x, [&](IdxX idx_x) {
        double const x = ddc::coordinate(idx_x);
        double const y_min(ymin);
        double const y_max(ymax);
        derivs_ymin_extracted(idx_x) = 2 * y_min + 2 * x;
        derivs_ymax_extracted(idx_x) = 2 * y_max + 2 * x;
    });

    ddc::for_each(deriv_block_xy, [&](Idx<dX, dY> idx_dxdy) {
        Idx<dX> idx_dx(idx_dxdy);
        Idx<dY> idx_dy(idx_dxdy);
        for (IdxX idx_x : idx_range_slice_dx) {
            for (IdxY idx_y : idx_range_slice_dy) {
                derivs_xy_extracted(
                        idx_dx - deriv_block_x.front(),
                        idx_dy - deriv_block_y.front(),
                        idx_range_slice_dx.distance_from_front(idx_x).value(),
                        idx_range_slice_dy.distance_from_front(idx_y).value())
                        = 2;
            }
        }
    });

    // Instantiate a spline representation builder -----------------------------------------------
    ddc::SplineBuilder2D<
            Kokkos::DefaultHostExecutionSpace,
            typename Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesX,
            BSplinesY,
            GridX,
            GridY,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::SplineSolver::LAPACK>
            builder(idx_range_xy);


    host_t<DFieldMem<IdxRange<BSplinesX, BSplinesY>>> function_coef_alloc(
            builder.batched_spline_domain(idx_range_xy));
    host_t<DField<IdxRange<BSplinesX, BSplinesY>>> function_coef = get_field(function_coef_alloc);

    // Copy the function values and its derivatives in the correct types -------------------------
    /*
        The SplineBuilder2D works with different types for the derivatives. 
        It words better for layout_right instead of layout_stride. 
    */
    Idx<DerivX> first_deriv_x(1);
    IdxStep<DerivX> n_deriv_x(1);
    IdxRange<DerivX> idx_range_deriv_x(first_deriv_x, n_deriv_x);

    Idx<DerivY> first_deriv_y(1);
    IdxStep<DerivY> n_deriv_y(1);
    IdxRange<DerivY> idx_range_deriv_y(first_deriv_y, n_deriv_y);

    Idx<DerivX, DerivY> first_deriv_x_deriv_y(first_deriv_x, first_deriv_y);
    IdxRange<DerivX, DerivY> idx_range_deriv_x_deriv_y(idx_range_deriv_x, idx_range_deriv_y);

    IdxRange<DerivX, GridY> idx_range_deriv_x_y(idx_range_deriv_x, idx_range_y);
    IdxRange<GridX, DerivY> idx_range_x_deriv_y(idx_range_x, idx_range_deriv_y);

    // --- Allocate memory
    host_t<DFieldMem<IdxRange<GridX, GridY>>> function_alloc(idx_range_xy);

    host_t<DFieldMem<IdxRange<DerivX, GridY>>> derivs_xmin_alloc(idx_range_deriv_x_y);
    host_t<DFieldMem<IdxRange<DerivX, GridY>>> derivs_xmax_alloc(idx_range_deriv_x_y);
    host_t<DFieldMem<IdxRange<GridX, DerivY>>> derivs_ymin_alloc(idx_range_x_deriv_y);
    host_t<DFieldMem<IdxRange<GridX, DerivY>>> derivs_ymax_alloc(idx_range_x_deriv_y);

    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> derivs_xy_min_min_alloc(idx_range_deriv_x_deriv_y);
    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> derivs_xy_max_min_alloc(idx_range_deriv_x_deriv_y);
    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> derivs_xy_min_max_alloc(idx_range_deriv_x_deriv_y);
    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> derivs_xy_max_max_alloc(idx_range_deriv_x_deriv_y);

    // --- Create Fields with the correct type
    host_t<DField<IdxRange<GridX, GridY>>> function(function_alloc);

    host_t<DField<IdxRange<DerivX, GridY>>> derivs_xmin(derivs_xmin_alloc);
    host_t<DField<IdxRange<DerivX, GridY>>> derivs_xmax(derivs_xmax_alloc);
    host_t<DField<IdxRange<GridX, DerivY>>> derivs_ymin(derivs_ymin_alloc);
    host_t<DField<IdxRange<GridX, DerivY>>> derivs_ymax(derivs_ymax_alloc);

    host_t<DField<IdxRange<DerivX, DerivY>>> derivs_xy_min_min(derivs_xy_min_min_alloc);
    host_t<DField<IdxRange<DerivX, DerivY>>> derivs_xy_max_min(derivs_xy_max_min_alloc);
    host_t<DField<IdxRange<DerivX, DerivY>>> derivs_xy_min_max(derivs_xy_min_max_alloc);
    host_t<DField<IdxRange<DerivX, DerivY>>> derivs_xy_max_max(derivs_xy_max_max_alloc);

    // --- Copy the values from the DerivField
    ddc::for_each(idx_range_xy, [&](IdxXY idx_xy) {
        function(idx_xy) = function_extracted(idx_xy);
    });

    ddc::for_each(idx_range_deriv_x_y, [&](Idx<DerivX, GridY> idx_deriv_x_y) {
        derivs_xmin(idx_deriv_x_y) = derivs_xmin_extracted(IdxY(idx_deriv_x_y));
        derivs_xmax(idx_deriv_x_y) = derivs_xmax_extracted(IdxY(idx_deriv_x_y));
    });

    ddc::for_each(idx_range_x_deriv_y, [&](Idx<GridX, DerivY> idx_x_deriv_y) {
        derivs_ymin(idx_x_deriv_y) = derivs_ymin_extracted(IdxX(idx_x_deriv_y));
        derivs_ymax(idx_x_deriv_y) = derivs_ymax_extracted(IdxX(idx_x_deriv_y));
    });

    derivs_xy_min_min(first_deriv_x_deriv_y) = derivs_xy_extracted(
            deriv_block_x.front() - deriv_block_x.front(), // IdxStep<dX>(0)
            deriv_block_y.front() - deriv_block_y.front(), // IdxStep<dY>(0)
            idx_range_slice_dx.distance_from_front(idx_range_slice_dx.front()).value(), // int(0)
            idx_range_slice_dy.distance_from_front(idx_range_slice_dy.front()).value()); // int(0)

    derivs_xy_max_min(first_deriv_x_deriv_y) = derivs_xy_extracted(
            deriv_block_x.back() - deriv_block_x.front(),
            deriv_block_y.front() - deriv_block_y.front(),
            idx_range_slice_dx.distance_from_front(idx_range_slice_dx.back()).value(),
            idx_range_slice_dy.distance_from_front(idx_range_slice_dy.front()).value());

    derivs_xy_min_max(first_deriv_x_deriv_y) = derivs_xy_extracted(
            deriv_block_x.front() - deriv_block_x.front(),
            deriv_block_y.back() - deriv_block_y.front(),
            idx_range_slice_dx.distance_from_front(idx_range_slice_dx.front()).value(),
            idx_range_slice_dy.distance_from_front(idx_range_slice_dy.back()).value());

    derivs_xy_max_max(first_deriv_x_deriv_y) = derivs_xy_extracted(
            deriv_block_x.back() - deriv_block_x.front(), // IdxStep<dX>(1)
            deriv_block_y.back() - deriv_block_y.front(), // IdxStep<dY>(1)
            idx_range_slice_dx.distance_from_front(idx_range_slice_dx.back()).value(), // int(1)
            idx_range_slice_dy.distance_from_front(idx_range_slice_dy.back()).value()); // int(1)

    // --- Build the spline representation
    builder(function_coef,
            get_const_field(function),
            std::optional(get_const_field(derivs_xmin)),
            std::optional(get_const_field(derivs_xmax)),
            std::optional(get_const_field(derivs_ymin)),
            std::optional(get_const_field(derivs_ymax)),
            std::optional(get_const_field(derivs_xy_min_min)),
            std::optional(get_const_field(derivs_xy_max_min)),
            std::optional(get_const_field(derivs_xy_min_max)),
            std::optional(get_const_field(derivs_xy_max_max)));


    // Test the spline representation ------------------------------------------------------------
    /*
        We verify here that the built spline representation is the expected one. 
        The spline is exact on the interpolation points. Its boundary derivatives 
        and corner cross-derivatives are also exact. 
    */
    // --- Define an evaluator
    ddc::ConstantExtrapolationRule<X, Y> bc_xmin(xmin, ymin, ymax);
    ddc::ConstantExtrapolationRule<X, Y> bc_xmax(xmax, ymin, ymax);
    ddc::ConstantExtrapolationRule<Y, X> bc_ymin(ymin, xmin, xmax);
    ddc::ConstantExtrapolationRule<Y, X> bc_ymax(ymax, xmin, xmax);
    ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            typename Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesX,
            BSplinesY,
            GridX,
            GridY,
            ddc::ConstantExtrapolationRule<X, Y>,
            ddc::ConstantExtrapolationRule<X, Y>,
            ddc::ConstantExtrapolationRule<Y, X>,
            ddc::ConstantExtrapolationRule<Y, X>>
            evaluator(bc_xmin, bc_xmax, bc_ymin, bc_ymax);

    // --- Check the function values
    ddc::for_each(idx_range_xy, [&](IdxXY idx_xy) {
        Coord<X, Y> eval_coord = ddc::coordinate(idx_xy);
        double const x = ddc::select<X>(eval_coord);
        double const y = ddc::select<Y>(eval_coord);

        const double spline_value = evaluator(eval_coord, get_const_field(function_coef));
        const double expected_value = x * x + y * y + 2 * x * y;

        EXPECT_NEAR(spline_value, expected_value, 1e-14);
    });

    // --- Check the x-derivatives
    ddc::for_each(idx_range_deriv_x_y, [&](Idx<DerivX, GridY> idx_deriv_x_y) {
        double const y = ddc::coordinate(IdxY(idx_deriv_x_y));
        double const x_min(xmin);
        double const x_max(xmax);

        Coord<X, Y> eval_coord(x_min, y);
        double spline_value = evaluator.deriv_dim_1(eval_coord, get_const_field(function_coef));
        double expected_value = 2 * x_min + 2 * y;
        EXPECT_NEAR(spline_value, expected_value, 1e-14);

        eval_coord = Coord<X, Y>(x_max, y);
        spline_value = evaluator.deriv_dim_1(eval_coord, get_const_field(function_coef));
        expected_value = 2 * x_max + 2 * y;
        EXPECT_NEAR(spline_value, expected_value, 1e-14);
    });

    // --- Check the y-derivatives
    ddc::for_each(idx_range_x_deriv_y, [&](Idx<GridX, DerivY> idx_x_deriv_y) {
        double const x = ddc::coordinate(IdxX(idx_x_deriv_y));
        double const y_min(ymin);
        double const y_max(ymax);

        Coord<X, Y> eval_coord(x, y_min);
        double spline_value = evaluator.deriv_dim_2(eval_coord, get_const_field(function_coef));
        double expected_value = 2 * y_min + 2 * x;
        EXPECT_NEAR(spline_value, expected_value, 1e-14);

        eval_coord = Coord<X, Y>(x, y_max);
        spline_value = evaluator.deriv_dim_2(eval_coord, get_const_field(function_coef));
        expected_value = 2 * y_max + 2 * x;
        EXPECT_NEAR(spline_value, expected_value, 1e-14);
    });

    // --- Check the cross-derivatives
    double expected_value = 2;

    Coord<X, Y> eval_coord(xmin, ymin);
    double spline_value = evaluator.deriv_1_and_2(eval_coord, get_const_field(function_coef));
    EXPECT_NEAR(spline_value, expected_value, 5e-14);

    eval_coord = Coord<X, Y>(xmax, ymin);
    spline_value = evaluator.deriv_1_and_2(eval_coord, get_const_field(function_coef));
    EXPECT_NEAR(spline_value, expected_value, 5e-14);

    eval_coord = Coord<X, Y>(xmin, ymax);
    spline_value = evaluator.deriv_1_and_2(eval_coord, get_const_field(function_coef));
    EXPECT_NEAR(spline_value, expected_value, 5e-14);

    eval_coord = Coord<X, Y>(xmax, ymax);
    spline_value = evaluator.deriv_1_and_2(eval_coord, get_const_field(function_coef));
    EXPECT_NEAR(spline_value, expected_value, 5e-14);
}