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

using DerivX = ddc::Deriv<X>;
using DerivY = ddc::Deriv<Y>;

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

using ExecSpace = Kokkos::DefaultExecutionSpace;

void initialise_function(
        DerivField<double, IdxRange<DerivX, GridX, DerivY, GridY>> function_and_derivs,
        DField<IdxRange<GridX, GridY>> function)
{
    ddc::parallel_for_each(
            ExecSpace(),
            get_idx_range(function),
            KOKKOS_LAMBDA(Idx<GridX, GridY> const idx) {
                double const x = ddc::coordinate(Idx<GridX>(idx));
                double const y = ddc::coordinate(Idx<GridY>(idx));
                function(idx) = Kokkos::cos(2. / 3 * M_PI * x + 0.25) * Kokkos::sin(y);
                function_and_derivs(idx) = function(idx);
            });
}

void initialise_derivatives(
        DerivField<double, IdxRange<DerivX, GridX, DerivY, GridY>> function_and_derivs,
        DField<IdxRange<DerivX, GridY>> derivs_xmin,
        DField<IdxRange<DerivX, GridY>> derivs_xmax,
        DField<IdxRange<GridX, DerivY>> derivs_ymin,
        DField<IdxRange<GridX, DerivY>> derivs_ymax,
        DField<IdxRange<DerivX, DerivY>> derivs_xy_min_min,
        DField<IdxRange<DerivX, DerivY>> derivs_xy_max_min,
        DField<IdxRange<DerivX, DerivY>> derivs_xy_min_max,
        DField<IdxRange<DerivX, DerivY>> derivs_xy_max_max)
{
    IdxRange<GridX> idx_range_x(get_idx_range(derivs_ymin));
    IdxRange<GridY> idx_range_y(get_idx_range(derivs_xmin));

    Coord<X> x_min(ddc::coordinate(idx_range_x.front()));
    Coord<X> x_max(ddc::coordinate(idx_range_x.back()));
    Coord<Y> y_min(ddc::coordinate(idx_range_y.front()));
    Coord<Y> y_max(ddc::coordinate(idx_range_y.back()));

    Idx<DerivX> first_dx(1);
    Idx<DerivY> first_dy(1);

    Idx<GridX> idx_slice_xmin(idx_range_x.front());
    Idx<GridX> idx_slice_xmax(idx_range_x.back());
    Idx<GridY> idx_slice_ymin(idx_range_y.front());
    Idx<GridY> idx_slice_ymax(idx_range_y.back());

    Idx<DerivX, GridX> idx_deriv_xmin(first_dx, idx_slice_xmin);
    Idx<DerivX, GridX> idx_deriv_xmax(first_dx, idx_slice_xmax);
    Idx<DerivY, GridY> idx_deriv_ymin(first_dy, idx_slice_ymin);
    Idx<DerivY, GridY> idx_deriv_ymax(first_dy, idx_slice_ymax);

    Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_min_min(idx_deriv_xmin, idx_deriv_ymin);
    Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_max_min(idx_deriv_xmax, idx_deriv_ymin);
    Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_min_max(idx_deriv_xmin, idx_deriv_ymax);
    Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_max_max(idx_deriv_xmax, idx_deriv_ymax);

    ddc::parallel_for_each(
            ExecSpace(),
            idx_range_y,
            KOKKOS_LAMBDA(Idx<GridY> const idx_y) {
                double const y = ddc::coordinate(idx_y);
                derivs_xmin(first_dx, idx_y) = -2. / 3 * M_PI
                                               * Kokkos::sin(2. / 3 * M_PI * double(x_min) + 0.25)
                                               * Kokkos::sin(y);
                derivs_xmax(first_dx, idx_y) = -2. / 3 * M_PI
                                               * Kokkos::sin(2. / 3 * M_PI * double(x_max) + 0.25)
                                               * Kokkos::sin(y);
                function_and_derivs(idx_deriv_xmin, idx_y) = derivs_xmin(first_dx, idx_y);
                function_and_derivs(idx_deriv_xmax, idx_y) = derivs_xmax(first_dx, idx_y);
            });

    ddc::parallel_for_each(
            ExecSpace(),
            idx_range_x,
            KOKKOS_LAMBDA(Idx<GridX> const idx_x) {
                double const x = ddc::coordinate(idx_x);
                derivs_ymin(idx_x, first_dy)
                        = Kokkos::cos(2. / 3 * M_PI * x + 0.25) * Kokkos::cos(double(y_min));
                derivs_ymax(idx_x, first_dy)
                        = Kokkos::cos(2. / 3 * M_PI * x + 0.25) * Kokkos::cos(double(y_max));
                function_and_derivs(idx_deriv_ymin, idx_x) = derivs_ymin(idx_x, first_dy);
                function_and_derivs(idx_deriv_ymax, idx_x) = derivs_ymax(idx_x, first_dy);
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
}

void initialise_derivatives_hybrid_case(
        DerivField<double, IdxRange<DerivX, GridX, DerivY, GridY>> function_and_derivs,
        DField<IdxRange<DerivX, GridY>> derivs_xmin,
        DField<IdxRange<GridX, DerivY>> derivs_ymax,
        DField<IdxRange<DerivX, DerivY>> derivs_xy_min_max)
{
    IdxRange<GridX> idx_range_x(get_idx_range(derivs_ymax));
    IdxRange<GridY> idx_range_y(get_idx_range(derivs_xmin));

    Coord<X> x_min(ddc::coordinate(idx_range_x.front()));
    Coord<Y> y_max(ddc::coordinate(idx_range_y.back()));

    Idx<DerivX> first_dx(1);
    Idx<DerivY> first_dy(1);

    Idx<GridX> idx_slice_xmin(idx_range_x.front());
    Idx<GridY> idx_slice_ymax(idx_range_y.back());

    Idx<DerivX, GridX> idx_deriv_xmin(first_dx, idx_slice_xmin);
    Idx<DerivY, GridY> idx_deriv_ymax(first_dy, idx_slice_ymax);
    Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_min_max(idx_deriv_xmin, idx_deriv_ymax);

    ddc::parallel_for_each(
            ExecSpace(),
            idx_range_y,
            KOKKOS_LAMBDA(Idx<GridY> const idx_y) {
                double const y = ddc::coordinate(idx_y);
                derivs_xmin(first_dx, idx_y) = -2. / 3 * M_PI
                                               * Kokkos::sin(2. / 3 * M_PI * double(x_min) + 0.25)
                                               * Kokkos::sin(y);
                function_and_derivs(idx_deriv_xmin, idx_y) = derivs_xmin(first_dx, idx_y);
            });

    ddc::parallel_for_each(
            ExecSpace(),
            idx_range_x,
            KOKKOS_LAMBDA(Idx<GridX> const idx_x) {
                double const x = ddc::coordinate(idx_x);
                derivs_ymax(idx_x, first_dy)
                        = Kokkos::cos(2. / 3 * M_PI * x + 0.25) * Kokkos::cos(double(y_max));
                function_and_derivs(idx_deriv_ymax, idx_x) = derivs_ymax(idx_x, first_dy);
            });


    Kokkos::parallel_for(
            "init_cross-derivs",
            Kokkos::RangePolicy<ExecSpace>(0, 1),
            KOKKOS_LAMBDA(const int) {
                derivs_xy_min_max(first_dx, first_dy) = -2. / 3 * M_PI
                                                        * Kokkos::sin(2. / 3 * M_PI * x_min + 0.25)
                                                        * Kokkos::sin(y_max);
                function_and_derivs(idx_cross_deriv_min_max)
                        = derivs_xy_min_max(first_dx, first_dy);
            });
}
} // namespace



TEST(SplineBuilderDerivField2DTest, DDCBoundCondHermiteTest)
{
    using SplineInterpPointsX = ddcHelper::NonUniformInterpolationPoints<
            BSplinesX,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE>;
    using SplineInterpPointsY = ddcHelper::NonUniformInterpolationPoints<
            BSplinesY,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE>;

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
    DerivFieldMem<double, IdxRange<DerivX, GridX, DerivY, GridY>, 1>
            function_and_derivs_alloc(idx_range_xy, idx_range_slice_dx, idx_range_slice_dy);
    DerivField<double, IdxRange<DerivX, GridX, DerivY, GridY>> function_and_derivs(
            function_and_derivs_alloc);

    // --- Fields
    DFieldMem<IdxRange<GridX, GridY>> function_alloc(idx_range_xy);
    DField<IdxRange<GridX, GridY>> function(function_alloc);

    Idx<DerivX> first_dx(1);
    IdxRange<DerivX> idx_range_deriv_x(first_dx, IdxStep<DerivX>(1));

    Idx<DerivY> first_dy(1);
    IdxRange<DerivY> idx_range_deriv_y(first_dy, IdxStep<DerivY>(1));

    IdxRange<DerivX, GridY> idx_range_dx_y(idx_range_deriv_x, idx_range_y);
    IdxRange<GridX, DerivY> idx_range_x_dy(idx_range_x, idx_range_deriv_y);
    IdxRange<DerivX, DerivY> idx_range_dx_dy(idx_range_deriv_x, idx_range_deriv_y);

    DFieldMem<IdxRange<DerivX, GridY>> derivs_xmin_alloc(idx_range_dx_y);
    DFieldMem<IdxRange<DerivX, GridY>> derivs_xmax_alloc(idx_range_dx_y);
    DFieldMem<IdxRange<GridX, DerivY>> derivs_ymin_alloc(idx_range_x_dy);
    DFieldMem<IdxRange<GridX, DerivY>> derivs_ymax_alloc(idx_range_x_dy);

    DFieldMem<IdxRange<DerivX, DerivY>> derivs_xy_min_min_alloc(idx_range_dx_dy);
    DFieldMem<IdxRange<DerivX, DerivY>> derivs_xy_max_min_alloc(idx_range_dx_dy);
    DFieldMem<IdxRange<DerivX, DerivY>> derivs_xy_min_max_alloc(idx_range_dx_dy);
    DFieldMem<IdxRange<DerivX, DerivY>> derivs_xy_max_max_alloc(idx_range_dx_dy);

    DField<IdxRange<DerivX, GridY>> derivs_xmin = get_field(derivs_xmin_alloc);
    DField<IdxRange<DerivX, GridY>> derivs_xmax = get_field(derivs_xmax_alloc);
    DField<IdxRange<GridX, DerivY>> derivs_ymin = get_field(derivs_ymin_alloc);
    DField<IdxRange<GridX, DerivY>> derivs_ymax = get_field(derivs_ymax_alloc);

    DField<IdxRange<DerivX, DerivY>> derivs_xy_min_min(derivs_xy_min_min_alloc);
    DField<IdxRange<DerivX, DerivY>> derivs_xy_max_min(derivs_xy_max_min_alloc);
    DField<IdxRange<DerivX, DerivY>> derivs_xy_min_max(derivs_xy_min_max_alloc);
    DField<IdxRange<DerivX, DerivY>> derivs_xy_max_max(derivs_xy_max_max_alloc);

    // Initialise data ---------------------------------------------------------------------------
    initialise_function(function_and_derivs, function);
    initialise_derivatives(
            function_and_derivs,
            derivs_xmin,
            derivs_xmax,
            derivs_ymin,
            derivs_ymax,
            derivs_xy_min_min,
            derivs_xy_max_min,
            derivs_xy_min_max,
            derivs_xy_max_max);

    // Instantiate the spline builders -----------------------------------------------------------
    ddc::SplineBuilder2D<
            ExecSpace,
            typename ExecSpace::memory_space,
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

    SplineBuliderDerivField2D<
            ExecSpace,
            BSplinesX,
            BSplinesY,
            GridX,
            GridY,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE>
            apply_builder(builder);

    // Instantiate splines -----------------------------------------------------------------------
    IdxRange<BSplinesX, BSplinesY> idx_range_spline = builder.batched_spline_domain(idx_range_xy);
    DFieldMem<IdxRange<BSplinesX, BSplinesY>> spline_deriv_field_alloc(idx_range_spline);
    DFieldMem<IdxRange<BSplinesX, BSplinesY>> spline_fields_alloc(idx_range_spline);

    DField<IdxRange<BSplinesX, BSplinesY>> spline_deriv_field(spline_deriv_field_alloc);
    DField<IdxRange<BSplinesX, BSplinesY>> spline_fields(spline_fields_alloc);

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

    ddc::host_for_each(idx_range_spline, [&](Idx<BSplinesX, BSplinesY> const& idx) {
        EXPECT_EQ(spline_deriv_field_host(idx), spline_fields_host(idx));
    });
}


TEST(SplineBuilderDerivField2DTest, DDCBoundCondGrevilleTest)
{
    using SplineInterpPointsX = ddc::GrevilleInterpolationPoints<
            BSplinesX,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE>;
    using SplineInterpPointsY = ddc::GrevilleInterpolationPoints<
            BSplinesY,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE>;

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

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridY>(SplineInterpPointsY::get_sampling<GridY>());

    IdxRange<GridX> idx_range_x(SplineInterpPointsX::template get_domain<GridX>());
    IdxRange<GridY> idx_range_y(SplineInterpPointsY::template get_domain<GridY>());
    IdxRange<GridX, GridY> idx_range_xy(idx_range_x, idx_range_y);

    IdxRangeSlice<GridX>
            idx_range_slice_dx(idx_range_x.front(), IdxStep<GridX>(0), idx_range_x.extents() - 1);
    IdxRangeSlice<GridY>
            idx_range_slice_dy(idx_range_y.front(), IdxStep<GridY>(0), idx_range_y.extents() - 1);

    // Instantiate data --------------------------------------------------------------------------
    // --- DerivField
    DerivFieldMem<double, IdxRange<DerivX, GridX, DerivY, GridY>, 0>
            function_and_derivs_alloc(idx_range_xy, idx_range_slice_dx, idx_range_slice_dy);
    DerivField<double, IdxRange<DerivX, GridX, DerivY, GridY>> function_and_derivs(
            function_and_derivs_alloc);

    // --- Fields
    DFieldMem<IdxRange<GridX, GridY>> function_alloc(idx_range_xy);
    DField<IdxRange<GridX, GridY>> function(function_alloc);

    // Initialise data ---------------------------------------------------------------------------
    initialise_function(function_and_derivs, function);

    // Instantiate the spline builders -----------------------------------------------------------
    ddc::SplineBuilder2D<
            ExecSpace,
            typename ExecSpace::memory_space,
            BSplinesX,
            BSplinesY,
            GridX,
            GridY,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::SplineSolver::LAPACK>
            builder(idx_range_xy);

    SplineBuliderDerivField2D<
            ExecSpace,
            BSplinesX,
            BSplinesY,
            GridX,
            GridY,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE>
            apply_builder(builder);

    // Instantiate splines -----------------------------------------------------------------------
    IdxRange<BSplinesX, BSplinesY> idx_range_spline = builder.batched_spline_domain(idx_range_xy);
    DFieldMem<IdxRange<BSplinesX, BSplinesY>> spline_deriv_field_alloc(idx_range_spline);
    DFieldMem<IdxRange<BSplinesX, BSplinesY>> spline_fields_alloc(idx_range_spline);

    DField<IdxRange<BSplinesX, BSplinesY>> spline_deriv_field(spline_deriv_field_alloc);
    DField<IdxRange<BSplinesX, BSplinesY>> spline_fields(spline_fields_alloc);

    // Build splines and compare -----------------------------------------------------------------
    apply_builder(spline_deriv_field, function_and_derivs);
    builder(spline_fields, get_const_field(function));

    auto spline_deriv_field_host = ddc::create_mirror_view_and_copy(spline_deriv_field);
    auto spline_fields_host = ddc::create_mirror_view_and_copy(spline_fields);

    ddc::host_for_each(idx_range_spline, [&](Idx<BSplinesX, BSplinesY> const& idx) {
        EXPECT_EQ(spline_deriv_field_host(idx), spline_fields_host(idx));
    });
}



TEST(SplineBuilderDerivField2DTest, HybridDDCBoundCondTest)
{
    using SplineInterpPointsX = ddcHelper::NonUniformInterpolationPoints<
            BSplinesX,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::GREVILLE>;
    using SplineInterpPointsY = ddcHelper::NonUniformInterpolationPoints<
            BSplinesY,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::HERMITE>;

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

    std::vector<Coord<X>> interpolation_points_x;
    for (int i(0); i < x_ncells.value(); i++) {
        interpolation_points_x.push_back(break_points_x[i]);
    }
    interpolation_points_x.push_back(
            1. / 3 * (break_points_x[x_ncells.value() - 1] + 2 * break_points_x[x_ncells.value()]));
    interpolation_points_x.push_back(break_points_x[x_ncells.value()]);

    std::vector<Coord<Y>> interpolation_points_y;
    interpolation_points_y.push_back(break_points_y[0]);
    interpolation_points_y.push_back(1. / 3 * (break_points_y[1] + 2 * break_points_y[0]));
    for (int i(1); i < y_ncells.value() + 1; i++) {
        interpolation_points_y.push_back(break_points_y[i]);
    }

    ddc::init_discrete_space<BSplinesX>(break_points_x);
    ddc::init_discrete_space<BSplinesY>(break_points_y);

    ddc::init_discrete_space<GridX>(interpolation_points_x);
    ddc::init_discrete_space<GridY>(interpolation_points_y);

    IdxRange<GridX> idx_range_x(SplineInterpPointsX::template get_domain<GridX>());
    IdxRange<GridY> idx_range_y(SplineInterpPointsY::template get_domain<GridY>());
    IdxRange<GridX, GridY> idx_range_xy(idx_range_x, idx_range_y);

    IdxRangeSlice<GridX>
            idx_range_slice_dx(idx_range_x.front(), IdxStep<GridX>(1), idx_range_x.extents() - 1);
    IdxRangeSlice<GridY>
            idx_range_slice_dy(idx_range_y.back(), IdxStep<GridY>(1), idx_range_y.extents() - 1);

    // Instantiate data --------------------------------------------------------------------------
    // --- DerivField
    DerivFieldMem<double, IdxRange<DerivX, GridX, DerivY, GridY>, 1>
            function_and_derivs_alloc(idx_range_xy, idx_range_slice_dx, idx_range_slice_dy);
    DerivField<double, IdxRange<DerivX, GridX, DerivY, GridY>> function_and_derivs(
            function_and_derivs_alloc);

    // --- Fields
    DFieldMem<IdxRange<GridX, GridY>> function_alloc(idx_range_xy);
    DField<IdxRange<GridX, GridY>> function(function_alloc);

    Idx<DerivX> first_dx(1);
    IdxRange<DerivX> idx_range_deriv_x(first_dx, IdxStep<DerivX>(1));

    Idx<DerivY> first_dy(1);
    IdxRange<DerivY> idx_range_deriv_y(first_dy, IdxStep<DerivY>(1));

    IdxRange<DerivX, GridY> idx_range_dx_y(idx_range_deriv_x, idx_range_y);
    IdxRange<GridX, DerivY> idx_range_x_dy(idx_range_x, idx_range_deriv_y);
    IdxRange<DerivX, DerivY> idx_range_dx_dy(idx_range_deriv_x, idx_range_deriv_y);

    DFieldMem<IdxRange<DerivX, GridY>> derivs_xmin_alloc(idx_range_dx_y);
    DFieldMem<IdxRange<GridX, DerivY>> derivs_ymax_alloc(idx_range_x_dy);
    DFieldMem<IdxRange<DerivX, DerivY>> derivs_xy_min_max_alloc(idx_range_dx_dy);

    DField<IdxRange<DerivX, GridY>> derivs_xmin = get_field(derivs_xmin_alloc);
    DField<IdxRange<GridX, DerivY>> derivs_ymax = get_field(derivs_ymax_alloc);
    DField<IdxRange<DerivX, DerivY>> derivs_xy_min_max(derivs_xy_min_max_alloc);

    // Initialise data ---------------------------------------------------------------------------
    initialise_function(function_and_derivs, function);
    initialise_derivatives_hybrid_case(
            function_and_derivs,
            derivs_xmin,
            derivs_ymax,
            derivs_xy_min_max);

    // Instantiate the spline builders -----------------------------------------------------------
    ddc::SplineBuilder2D<
            ExecSpace,
            typename ExecSpace::memory_space,
            BSplinesX,
            BSplinesY,
            GridX,
            GridY,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::HERMITE,
            ddc::SplineSolver::LAPACK>
            builder(idx_range_xy);

    SplineBuliderDerivField2D<
            ExecSpace,
            BSplinesX,
            BSplinesY,
            GridX,
            GridY,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::HERMITE>
            apply_builder(builder);

    // Instantiate splines -----------------------------------------------------------------------
    IdxRange<BSplinesX, BSplinesY> idx_range_spline = builder.batched_spline_domain(idx_range_xy);
    DFieldMem<IdxRange<BSplinesX, BSplinesY>> spline_deriv_field_alloc(idx_range_spline);
    DFieldMem<IdxRange<BSplinesX, BSplinesY>> spline_fields_alloc(idx_range_spline);

    DField<IdxRange<BSplinesX, BSplinesY>> spline_deriv_field(spline_deriv_field_alloc);
    DField<IdxRange<BSplinesX, BSplinesY>> spline_fields(spline_fields_alloc);

    // Build splines and compare -----------------------------------------------------------------
    apply_builder(spline_deriv_field, function_and_derivs);

    builder(spline_fields,
            get_const_field(function),
            std::optional(get_const_field(derivs_xmin)),
            std::optional<DConstField<IdxRange<DerivX, GridY>>> {std::nullopt},
            std::optional<DConstField<IdxRange<GridX, DerivY>>> {std::nullopt},
            std::optional(get_const_field(derivs_ymax)),
            std::optional<DConstField<IdxRange<DerivX, DerivY>>> {std::nullopt},
            std::optional<DConstField<IdxRange<DerivX, DerivY>>> {std::nullopt},
            std::optional(get_const_field(derivs_xy_min_max)),
            std::optional<DConstField<IdxRange<DerivX, DerivY>>> {std::nullopt});

    auto spline_deriv_field_host = ddc::create_mirror_view_and_copy(spline_deriv_field);
    auto spline_fields_host = ddc::create_mirror_view_and_copy(spline_fields);

    ddc::host_for_each(idx_range_spline, [&](Idx<BSplinesX, BSplinesY> const& idx) {
        EXPECT_EQ(spline_deriv_field_host(idx), spline_fields_host(idx));
    });
}