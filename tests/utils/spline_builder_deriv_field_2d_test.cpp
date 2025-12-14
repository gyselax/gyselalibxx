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

// Grids for the ddc::BoundCond::HERMITE case.
struct GridXH : NonUniformGridBase<X>
{
};
struct GridYH : NonUniformGridBase<Y>
{
};

struct BSplinesXH : ddc::NonUniformBSplines<X, 3>
{
};
struct BSplinesYH : ddc::NonUniformBSplines<Y, 3>
{
};

// Grids for the ddc::BoundCond::GREVILLE case.
struct GridXG : NonUniformGridBase<X>
{
};
struct GridYG : NonUniformGridBase<Y>
{
};

struct BSplinesXG : ddc::NonUniformBSplines<X, 3>
{
};
struct BSplinesYG : ddc::NonUniformBSplines<Y, 3>
{
};


using ExecSpace = Kokkos::DefaultExecutionSpace;

} // namespace



TEST(SplineBuilderDerivField2DTest, DDCBoundCondHermiteTest)
{
    using DerivX = ddc::Deriv<X>;
    using DerivY = ddc::Deriv<Y>;

    using SplineInterpPointsX = ddcHelper::NonUniformInterpolationPoints<
            BSplinesXH,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE>;
    using SplineInterpPointsY = ddcHelper::NonUniformInterpolationPoints<
            BSplinesYH,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE>;


    // Initialise the mesh -----------------------------------------------------------------------
    const Coord<X> x_min(0.0);
    const Coord<X> x_max(1.0);
    const IdxStep<GridXH> x_ncells(8);

    const Coord<Y> y_min(0.0);
    const Coord<Y> y_max(1.0);
    const IdxStep<GridYH> y_ncells(8);

    std::vector<Coord<X>> break_points_x
            = build_random_non_uniform_break_points(x_min, x_max, x_ncells);
    std::vector<Coord<Y>> break_points_y
            = build_random_non_uniform_break_points(y_min, y_max, y_ncells);

    ddc::init_discrete_space<BSplinesXH>(break_points_x);
    ddc::init_discrete_space<BSplinesYH>(break_points_y);

    ddc::init_discrete_space<GridXH>(break_points_x);
    ddc::init_discrete_space<GridYH>(break_points_y);

    IdxRange<GridXH> idx_range_x(SplineInterpPointsX::template get_domain<GridXH>());
    IdxRange<GridYH> idx_range_y(SplineInterpPointsY::template get_domain<GridYH>());
    IdxRange<GridXH, GridYH> idx_range_xy(idx_range_x, idx_range_y);

    IdxRangeSlice<GridXH>
            idx_range_slice_dx(idx_range_x.front(), IdxStep<GridXH>(2), idx_range_x.extents() - 1);
    IdxRangeSlice<GridYH>
            idx_range_slice_dy(idx_range_y.front(), IdxStep<GridYH>(2), idx_range_y.extents() - 1);

    // Instantiate data --------------------------------------------------------------------------
    // --- DerivField
    DerivFieldMem<double, IdxRange<DerivX, GridXH, DerivY, GridYH>, 1>
            function_and_derivs_alloc(idx_range_xy, idx_range_slice_dx, idx_range_slice_dy);
    DerivField<double, IdxRange<DerivX, GridXH, DerivY, GridYH>> function_and_derivs(
            function_and_derivs_alloc);

    // --- Fields
    DFieldMem<IdxRange<GridXH, GridYH>> function_alloc(idx_range_xy);
    DField<IdxRange<GridXH, GridYH>> function(function_alloc);

    Idx<DerivX> first_dx(1);
    IdxRange<DerivX> idx_range_deriv_x(first_dx, IdxStep<DerivX>(1));

    Idx<DerivY> first_dy(1);
    IdxRange<DerivY> idx_range_deriv_y(first_dy, IdxStep<DerivY>(1));

    IdxRange<DerivX, GridYH> idx_range_dx_y(idx_range_deriv_x, idx_range_y);
    IdxRange<GridXH, DerivY> idx_range_x_dy(idx_range_x, idx_range_deriv_y);
    IdxRange<DerivX, DerivY> idx_range_dx_dy(idx_range_deriv_x, idx_range_deriv_y);

    DFieldMem<IdxRange<DerivX, GridYH>> derivs_xmin_alloc(idx_range_dx_y);
    DFieldMem<IdxRange<DerivX, GridYH>> derivs_xmax_alloc(idx_range_dx_y);
    DFieldMem<IdxRange<GridXH, DerivY>> derivs_ymin_alloc(idx_range_x_dy);
    DFieldMem<IdxRange<GridXH, DerivY>> derivs_ymax_alloc(idx_range_x_dy);

    DFieldMem<IdxRange<DerivX, DerivY>> derivs_xy_min_min_alloc(idx_range_dx_dy);
    DFieldMem<IdxRange<DerivX, DerivY>> derivs_xy_max_min_alloc(idx_range_dx_dy);
    DFieldMem<IdxRange<DerivX, DerivY>> derivs_xy_min_max_alloc(idx_range_dx_dy);
    DFieldMem<IdxRange<DerivX, DerivY>> derivs_xy_max_max_alloc(idx_range_dx_dy);

    DField<IdxRange<DerivX, GridYH>> derivs_xmin = get_field(derivs_xmin_alloc);
    DField<IdxRange<DerivX, GridYH>> derivs_xmax = get_field(derivs_xmax_alloc);
    DField<IdxRange<GridXH, DerivY>> derivs_ymin = get_field(derivs_ymin_alloc);
    DField<IdxRange<GridXH, DerivY>> derivs_ymax = get_field(derivs_ymax_alloc);

    DField<IdxRange<DerivX, DerivY>> derivs_xy_min_min(derivs_xy_min_min_alloc);
    DField<IdxRange<DerivX, DerivY>> derivs_xy_max_min(derivs_xy_max_min_alloc);
    DField<IdxRange<DerivX, DerivY>> derivs_xy_min_max(derivs_xy_min_max_alloc);
    DField<IdxRange<DerivX, DerivY>> derivs_xy_max_max(derivs_xy_max_max_alloc);

    // Initialise data ---------------------------------------------------------------------------
    Idx<GridXH> idx_slice_xmin(idx_range_slice_dx.front());
    Idx<GridXH> idx_slice_xmax(idx_range_slice_dx.back());
    Idx<GridYH> idx_slice_ymin(idx_range_slice_dy.front());
    Idx<GridYH> idx_slice_ymax(idx_range_slice_dy.back());

    Idx<DerivX, GridXH> idx_deriv_xmin(first_dx, idx_slice_xmin);
    Idx<DerivX, GridXH> idx_deriv_xmax(first_dx, idx_slice_xmax);
    Idx<DerivY, GridYH> idx_deriv_ymin(first_dy, idx_slice_ymin);
    Idx<DerivY, GridYH> idx_deriv_ymax(first_dy, idx_slice_ymax);

    Idx<DerivX, GridXH, DerivY, GridYH> idx_cross_deriv_min_min(idx_deriv_xmin, idx_deriv_ymin);
    Idx<DerivX, GridXH, DerivY, GridYH> idx_cross_deriv_max_min(idx_deriv_xmax, idx_deriv_ymin);
    Idx<DerivX, GridXH, DerivY, GridYH> idx_cross_deriv_min_max(idx_deriv_xmin, idx_deriv_ymax);
    Idx<DerivX, GridXH, DerivY, GridYH> idx_cross_deriv_max_max(idx_deriv_xmax, idx_deriv_ymax);

    ddc::parallel_for_each(
            ExecSpace(),
            idx_range_xy,
            KOKKOS_LAMBDA(Idx<GridXH, GridYH> const idx) {
                double const x = ddc::coordinate(Idx<GridXH>(idx));
                double const y = ddc::coordinate(Idx<GridYH>(idx));
                function(idx) = Kokkos::cos(2. / 3 * M_PI * x + 0.25) * Kokkos::sin(y);
                function_and_derivs.get_values_field()(idx) = function(idx);
            });

    ddc::parallel_for_each(
            ExecSpace(),
            idx_range_y,
            KOKKOS_LAMBDA(Idx<GridYH> const idx_y) {
                double const y = ddc::coordinate(idx_y);
                derivs_xmin(first_dx, idx_y) = -2. / 3 * M_PI
                                               * Kokkos::sin(2. / 3 * M_PI * double(x_min) + 0.25)
                                               * Kokkos::sin(y);
                derivs_xmax(first_dx, idx_y) = -2. / 3 * M_PI
                                               * Kokkos::sin(2. / 3 * M_PI * double(x_max) + 0.25)
                                               * Kokkos::sin(y);
                function_and_derivs[idx_deriv_xmin](idx_y) = derivs_xmin(first_dx, idx_y);
                function_and_derivs[idx_deriv_xmax](idx_y) = derivs_xmax(first_dx, idx_y);
            });

    ddc::parallel_for_each(
            ExecSpace(),
            idx_range_x,
            KOKKOS_LAMBDA(Idx<GridXH> const idx_x) {
                double const x = ddc::coordinate(idx_x);
                derivs_ymin(idx_x, first_dy)
                        = Kokkos::cos(2. / 3 * M_PI * x + 0.25) * Kokkos::cos(double(y_min));
                derivs_ymax(idx_x, first_dy)
                        = Kokkos::cos(2. / 3 * M_PI * x + 0.25) * Kokkos::cos(double(y_max));
                function_and_derivs[idx_deriv_ymin](idx_x) = derivs_ymin(idx_x, first_dy);
                function_and_derivs[idx_deriv_ymax](idx_x) = derivs_ymax(idx_x, first_dy);
            });

    Kokkos::parallel_for(
            "init_cross-derivs",
            Kokkos::RangePolicy<ExecSpace>(0, 1),
            KOKKOS_LAMBDA(const int) {
                derivs_xy_min_min(first_dx, first_dy) = -2. / 3 * M_PI
                                                        * Kokkos::sin(2. / 3 * M_PI * x_min + 0.25)
                                                        * Kokkos::sin(y_min);
                derivs_xy_max_min(first_dx, first_dy) = -2. / 3 * M_PI
                                                        * Kokkos::sin(2. / 3 * M_PI * x_max + 0.25)
                                                        * Kokkos::sin(y_min);
                derivs_xy_min_max(first_dx, first_dy) = -2. / 3 * M_PI
                                                        * Kokkos::sin(2. / 3 * M_PI * x_min + 0.25)
                                                        * Kokkos::sin(y_max);
                derivs_xy_max_max(first_dx, first_dy) = -2. / 3 * M_PI
                                                        * Kokkos::sin(2. / 3 * M_PI * x_max + 0.25)
                                                        * Kokkos::sin(y_max);

                function_and_derivs(idx_cross_deriv_min_min)
                        = derivs_xy_min_min(first_dx, first_dy);
                function_and_derivs(idx_cross_deriv_max_min)
                        = derivs_xy_max_min(first_dx, first_dy);
                function_and_derivs(idx_cross_deriv_min_max)
                        = derivs_xy_min_max(first_dx, first_dy);
                function_and_derivs(idx_cross_deriv_max_max)
                        = derivs_xy_max_max(first_dx, first_dy);
            });

    // Instantiate the spline builders -----------------------------------------------------------
    ddc::SplineBuilder2D<
            ExecSpace,
            typename ExecSpace::memory_space,
            BSplinesXH,
            BSplinesYH,
            GridXH,
            GridYH,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::SplineSolver::LAPACK>
            builder(idx_range_xy);

    SplineBuliderDerivField2D<
            ExecSpace,
            BSplinesXH,
            BSplinesYH,
            GridXH,
            GridYH,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE>
            apply_builder(builder);

    // Instantiate splines -----------------------------------------------------------------------
    IdxRange<BSplinesXH, BSplinesYH> idx_range_spline = builder.batched_spline_domain(idx_range_xy);
    DFieldMem<IdxRange<BSplinesXH, BSplinesYH>> spline_deriv_field_alloc(idx_range_spline);
    DFieldMem<IdxRange<BSplinesXH, BSplinesYH>> spline_fields_alloc(idx_range_spline);

    DField<IdxRange<BSplinesXH, BSplinesYH>> spline_deriv_field(spline_deriv_field_alloc);
    DField<IdxRange<BSplinesXH, BSplinesYH>> spline_fields(spline_fields_alloc);

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

    auto spline_deriv_field_host = ddc::create_mirror_view_and_copy(spline_deriv_field);
    auto spline_fields_host = ddc::create_mirror_view_and_copy(spline_fields);

    ddc::for_each(idx_range_spline, [&](Idx<BSplinesXH, BSplinesYH> const& idx) {
        EXPECT_EQ(spline_deriv_field_host(idx), spline_fields_host(idx));
    });
}


TEST(SplineBuilderDerivField2DTest, DDCBoundCondGrevilleTest)
{
    using DerivX = ddc::Deriv<X>;
    using DerivY = ddc::Deriv<Y>;

    using SplineInterpPointsX = ddc::GrevilleInterpolationPoints<
            BSplinesXG,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE>;
    using SplineInterpPointsY = ddc::GrevilleInterpolationPoints<
            BSplinesYG,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE>;


    // Initialise the mesh -----------------------------------------------------------------------
    const Coord<X> x_min(0.0);
    const Coord<X> x_max(1.0);
    const IdxStep<GridXG> x_ncells(8);

    const Coord<Y> y_min(0.0);
    const Coord<Y> y_max(1.0);
    const IdxStep<GridYG> y_ncells(8);

    std::vector<Coord<X>> break_points_x
            = build_random_non_uniform_break_points(x_min, x_max, x_ncells);
    std::vector<Coord<Y>> break_points_y
            = build_random_non_uniform_break_points(y_min, y_max, y_ncells);

    ddc::init_discrete_space<BSplinesXG>(break_points_x);
    ddc::init_discrete_space<BSplinesYG>(break_points_y);

    ddc::init_discrete_space<GridXG>(SplineInterpPointsX::get_sampling<GridXG>());
    ddc::init_discrete_space<GridYG>(SplineInterpPointsY::get_sampling<GridYG>());

    IdxRange<GridXG> idx_range_x(SplineInterpPointsX::template get_domain<GridXG>());
    IdxRange<GridYG> idx_range_y(SplineInterpPointsY::template get_domain<GridYG>());
    IdxRange<GridXG, GridYG> idx_range_xy(idx_range_x, idx_range_y);

    IdxRangeSlice<GridXG>
            idx_range_slice_dx(idx_range_x.front(), IdxStep<GridXG>(0), idx_range_x.extents() - 1);
    IdxRangeSlice<GridYG>
            idx_range_slice_dy(idx_range_y.front(), IdxStep<GridYG>(0), idx_range_y.extents() - 1);

    // Instantiate data --------------------------------------------------------------------------
    // --- DerivField
    DerivFieldMem<double, IdxRange<DerivX, GridXG, DerivY, GridYG>, 0>
            function_and_derivs_alloc(idx_range_xy, idx_range_slice_dx, idx_range_slice_dy);
    DerivField<double, IdxRange<DerivX, GridXG, DerivY, GridYG>> function_and_derivs(
            function_and_derivs_alloc);

    // --- Fields
    DFieldMem<IdxRange<GridXG, GridYG>> function_alloc(idx_range_xy);
    DField<IdxRange<GridXG, GridYG>> function(function_alloc);

    // Initialise data ---------------------------------------------------------------------------
    ddc::parallel_for_each(
            ExecSpace(),
            idx_range_xy,
            KOKKOS_LAMBDA(Idx<GridXG, GridYG> const idx) {
                double const x = ddc::coordinate(Idx<GridXG>(idx));
                double const y = ddc::coordinate(Idx<GridYG>(idx));
                function(idx) = Kokkos::cos(2. / 3 * M_PI * x + 0.25) * Kokkos::sin(y);
                function_and_derivs.get_values_field()(idx) = function(idx);
            });

    // Instantiate the spline builders -----------------------------------------------------------
    ddc::SplineBuilder2D<
            ExecSpace,
            typename ExecSpace::memory_space,
            BSplinesXG,
            BSplinesYG,
            GridXG,
            GridYG,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::SplineSolver::LAPACK>
            builder(idx_range_xy);

    SplineBuliderDerivField2D<
            ExecSpace,
            BSplinesXG,
            BSplinesYG,
            GridXG,
            GridYG,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE>
            apply_builder(builder);

    // Instantiate splines -----------------------------------------------------------------------
    IdxRange<BSplinesXG, BSplinesYG> idx_range_spline = builder.batched_spline_domain(idx_range_xy);
    DFieldMem<IdxRange<BSplinesXG, BSplinesYG>> spline_deriv_field_alloc(idx_range_spline);
    DFieldMem<IdxRange<BSplinesXG, BSplinesYG>> spline_fields_alloc(idx_range_spline);

    DField<IdxRange<BSplinesXG, BSplinesYG>> spline_deriv_field(spline_deriv_field_alloc);
    DField<IdxRange<BSplinesXG, BSplinesYG>> spline_fields(spline_fields_alloc);

    // Build splines and compare -----------------------------------------------------------------
    apply_builder(spline_deriv_field, function_and_derivs);
    builder(spline_fields, get_const_field(function));

    auto spline_deriv_field_host = ddc::create_mirror_view_and_copy(spline_deriv_field);
    auto spline_fields_host = ddc::create_mirror_view_and_copy(spline_fields);

    ddc::for_each(idx_range_spline, [&](Idx<BSplinesXG, BSplinesYG> const& idx) {
        EXPECT_EQ(spline_deriv_field_host(idx), spline_fields_host(idx));
    });
}