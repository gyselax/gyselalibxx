// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "spline_builder_deriv_field_2d.hpp"

namespace {
struct X
{
    static bool constexpr PERIODIC = false;
};
struct Y
{
    static bool constexpr PERIODIC = false;
};

struct GridX : NonUniformGridBase<X>
{
};
struct GridY : NonUniformGridBase<Y>
{
};

struct BSplinesX : ddc::NonUniformBSplines<X, 3>
{
};
struct BSplinesY : ddc::NonUniformBSplines<Y, 3>
{
};

using DerivX = ddc::Deriv<X>;
using DerivY = ddc::Deriv<Y>;

using HostExecSpace = Kokkos::DefaultHostExecutionSpace;

using SplineInterpPointsX = ddcHelper::
        NonUniformInterpolationPoints<BSplinesX, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;
using SplineInterpPointsY = ddcHelper::
        NonUniformInterpolationPoints<BSplinesY, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;

using SplineBuilder = ddc::SplineBuilder2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesX,
        BSplinesY,
        GridX,
        GridY,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::SplineSolver::LAPACK>;

using SplineBuilderDerivField = SplineBuliderDerivField2D<
        HostExecSpace,
        BSplinesX,
        BSplinesY,
        GridX,
        GridY,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

} // namespace



TEST(SplineBuilderDerivField2DTest, CorrectSplineTest)
{
    // Initialise the mesh -----------------------------------------------------------------------
    const Coord<X> x_min(0.0);
    const Coord<X> x_max(1.0);
    const IdxStep<GridX> x_ncells(8);

    const Coord<Y> y_min(0.0);
    const Coord<Y> y_max(1.0);
    const IdxStep<GridY> y_ncells(8);

    std::vector<Coord<X>> break_points_x
            = build_random_non_uniform_break_points(x_min, x_max, x_ncells);
    std::vector<Coord<Y>> break_points_y
            = build_random_non_uniform_break_points(y_min, y_max, y_ncells);

    ddc::init_discrete_space<BSplinesX>(break_points_x);
    ddc::init_discrete_space<BSplinesY>(break_points_y);

    ddc::init_discrete_space<GridX>(break_points_x);
    ddc::init_discrete_space<GridY>(break_points_y);

    IdxRange<GridX> idx_range_x(SplineInterpPointsX::template get_domain<GridX>());
    IdxRange<GridY> idx_range_y(SplineInterpPointsY::template get_domain<GridY>());
    IdxRange<GridX, GridY> idx_range_xy(idx_range_x, idx_range_y);

    IdxRangeSlice<GridX>
            idx_range_slice_dx(idx_range_x.front(), IdxStep<GridX>(2), idx_range_x.extents() - 1);
    IdxRangeSlice<GridY>
            idx_range_slice_dy(idx_range_y.front(), IdxStep<GridY>(2), idx_range_y.extents() - 1);

    // Instantiate data --------------------------------------------------------------------------
    // --- DerivField
    host_t<DerivFieldMem<double, IdxRange<DerivX, GridX, DerivY, GridY>, 1>>
            function_and_derivs_alloc(idx_range_xy, idx_range_slice_dx, idx_range_slice_dy);
    host_t<DerivField<double, IdxRange<DerivX, GridX, DerivY, GridY>>> function_and_derivs(
            function_and_derivs_alloc);

    // --- Fields
    host_t<DFieldMem<IdxRange<GridX, GridY>>> function_alloc(idx_range_xy);
    host_t<DField<IdxRange<GridX, GridY>>> function(function_alloc);

    Idx<DerivX> first_dx(1);
    IdxRange<DerivX> idx_range_deriv_x(first_dx, IdxStep<DerivX>(1));

    Idx<DerivY> first_dy(1);
    IdxRange<DerivY> idx_range_deriv_y(first_dy, IdxStep<DerivY>(1));

    IdxRange<DerivX, GridY> idx_range_dx_y(idx_range_deriv_x, idx_range_y);
    IdxRange<GridX, DerivY> idx_range_x_dy(idx_range_x, idx_range_deriv_y);
    IdxRange<DerivX, DerivY> idx_range_dx_dy(idx_range_deriv_x, idx_range_deriv_y);

    host_t<DFieldMem<IdxRange<DerivX, GridY>>> derivs_xmin_alloc(idx_range_dx_y);
    host_t<DFieldMem<IdxRange<DerivX, GridY>>> derivs_xmax_alloc(idx_range_dx_y);
    host_t<DFieldMem<IdxRange<GridX, DerivY>>> derivs_ymin_alloc(idx_range_x_dy);
    host_t<DFieldMem<IdxRange<GridX, DerivY>>> derivs_ymax_alloc(idx_range_x_dy);

    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> derivs_xy_min_min_alloc(idx_range_dx_dy);
    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> derivs_xy_max_min_alloc(idx_range_dx_dy);
    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> derivs_xy_min_max_alloc(idx_range_dx_dy);
    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> derivs_xy_max_max_alloc(idx_range_dx_dy);

    host_t<DField<IdxRange<DerivX, GridY>>> derivs_xmin = get_field(derivs_xmin_alloc);
    host_t<DField<IdxRange<DerivX, GridY>>> derivs_xmax = get_field(derivs_xmax_alloc);
    host_t<DField<IdxRange<GridX, DerivY>>> derivs_ymin = get_field(derivs_ymin_alloc);
    host_t<DField<IdxRange<GridX, DerivY>>> derivs_ymax = get_field(derivs_ymax_alloc);

    host_t<DField<IdxRange<DerivX, DerivY>>> derivs_xy_min_min(derivs_xy_min_min_alloc);
    host_t<DField<IdxRange<DerivX, DerivY>>> derivs_xy_max_min(derivs_xy_max_min_alloc);
    host_t<DField<IdxRange<DerivX, DerivY>>> derivs_xy_min_max(derivs_xy_min_max_alloc);
    host_t<DField<IdxRange<DerivX, DerivY>>> derivs_xy_max_max(derivs_xy_max_max_alloc);

    // Initialise data ---------------------------------------------------------------------------
    Idx<GridX> idx_slice_xmin(idx_range_slice_dx.front());
    Idx<GridX> idx_slice_xmax(idx_range_slice_dx.back());
    Idx<GridY> idx_slice_ymin(idx_range_slice_dy.front());
    Idx<GridY> idx_slice_ymax(idx_range_slice_dy.back());

    Idx<DerivX, GridX> idx_deriv_xmin(first_dx, idx_slice_xmin);
    Idx<DerivX, GridX> idx_deriv_xmax(first_dx, idx_slice_xmax);
    Idx<DerivY, GridY> idx_deriv_ymin(first_dy, idx_slice_ymin);
    Idx<DerivY, GridY> idx_deriv_ymax(first_dy, idx_slice_ymax);

    Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_min_min(idx_deriv_xmin, idx_deriv_ymin);
    Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_max_min(idx_deriv_xmax, idx_deriv_ymin);
    Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_min_max(idx_deriv_xmin, idx_deriv_ymax);
    Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_max_max(idx_deriv_xmax, idx_deriv_ymax);

    ddc::parallel_for_each(
            HostExecSpace(),
            idx_range_xy,
            KOKKOS_LAMBDA(Idx<GridX, GridY> const idx) {
                double const x = ddc::coordinate(Idx<GridX>(idx));
                double const y = ddc::coordinate(Idx<GridY>(idx));
                function(idx) = std::cos(2. / 3 * M_PI * x + 0.25) * std::sin(y);
                function_and_derivs.get_values_field()(idx) = function(idx);
            });

    ddc::for_each(idx_range_y, [&](Idx<GridY> const idx_y) {
        double const y = ddc::coordinate(idx_y);
        derivs_xmin(first_dx, idx_y)
                = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * double(x_min) + 0.25) * std::sin(y);
        derivs_xmax(first_dx, idx_y)
                = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * double(x_max) + 0.25) * std::sin(y);
        function_and_derivs[idx_deriv_xmin](idx_y) = derivs_xmin(first_dx, idx_y);
        function_and_derivs[idx_deriv_xmax](idx_y) = derivs_xmax(first_dx, idx_y);
    });

    ddc::for_each(idx_range_x, [&](Idx<GridX> const idx_x) {
        double const x = ddc::coordinate(idx_x);
        derivs_ymin(idx_x, first_dy)
                = std::cos(2. / 3 * M_PI * x + 0.25) * std ::cos(double(y_min));
        derivs_ymax(idx_x, first_dy)
                = std::cos(2. / 3 * M_PI * x + 0.25) * std ::cos(double(y_max));
        function_and_derivs[idx_deriv_ymin](idx_x) = derivs_ymin(idx_x, first_dy);
        function_and_derivs[idx_deriv_ymax](idx_x) = derivs_ymax(idx_x, first_dy);
    });

    derivs_xy_min_min(first_dx, first_dy)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * x_min + 0.25) * std::sin(y_min);
    derivs_xy_max_min(first_dx, first_dy)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * x_max + 0.25) * std::sin(y_min);
    derivs_xy_min_max(first_dx, first_dy)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * x_min + 0.25) * std::sin(y_max);
    derivs_xy_max_max(first_dx, first_dy)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * x_max + 0.25) * std::sin(y_max);

    function_and_derivs(idx_cross_deriv_min_min) = derivs_xy_min_min(first_dx, first_dy);
    function_and_derivs(idx_cross_deriv_max_min) = derivs_xy_max_min(first_dx, first_dy);
    function_and_derivs(idx_cross_deriv_min_max) = derivs_xy_min_max(first_dx, first_dy);
    function_and_derivs(idx_cross_deriv_max_max) = derivs_xy_max_max(first_dx, first_dy);

    // Instantiate the spline builders -----------------------------------------------------------
    SplineBuilder builder(idx_range_xy);
    SplineBuilderDerivField apply_builder(builder);

    // Instantiate splines -----------------------------------------------------------------------
    IdxRange<BSplinesX, BSplinesY> idx_range_spline = builder.batched_spline_domain(idx_range_xy);
    host_t<DFieldMem<IdxRange<BSplinesX, BSplinesY>>> spline_deriv_field_alloc(idx_range_spline);
    host_t<DFieldMem<IdxRange<BSplinesX, BSplinesY>>> spline_fields_alloc(idx_range_spline);

    host_t<DField<IdxRange<BSplinesX, BSplinesY>>> spline_deriv_field(spline_deriv_field_alloc);
    host_t<DField<IdxRange<BSplinesX, BSplinesY>>> spline_fields(spline_fields_alloc);

    // Build splines and compare -----------------------------------------------------------------
    apply_builder(spline_deriv_field, function_and_derivs);

    builder(spline_fields,
            get_const_field(function),
            std::optional(get_const_field(derivs_xmin)),
            std::optional(get_const_field(derivs_xmax)),
            std::optional(get_const_field(derivs_ymin)),
            std::optional(get_const_field(derivs_ymax)),
            std::optional(get_const_field(derivs_xy_min_min)),
            std::optional(get_const_field(derivs_xy_max_min)),
            std::optional(get_const_field(derivs_xy_min_max)),
            std::optional(get_const_field(derivs_xy_max_max)));

    ddc::for_each(idx_range_spline, [&](Idx<BSplinesX, BSplinesY> const& idx) {
        EXPECT_EQ(spline_deriv_field(idx), spline_fields(idx));
    });
}