// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_non_periodic_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "mesh_builder.hpp"
#include "multipatch_field.hpp"
#include "multipatch_spline_builder_2d.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"


namespace {
using namespace non_periodic_non_uniform_2d_2patches;

template <int PatchIdx, ddc::BoundCond BC>
using SplineInterpPointsX = ddc::GrevilleInterpolationPoints<BSplinesX<PatchIdx>, BC, BC>;
template <int PatchIdx, ddc::BoundCond BC>
using SplineInterpPointsY = ddc::GrevilleInterpolationPoints<BSplinesY<PatchIdx>, BC, BC>;

// Operators
template <int PatchIdx, ddc::BoundCond BC>
using SplineXYBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX<PatchIdx>,
        BSplinesY<PatchIdx>,
        GridX<PatchIdx>,
        GridY<PatchIdx>,
        BC,
        BC,
        BC,
        BC,
        ddc::SplineSolver::LAPACK>;

template <ddc::BoundCond BC>
using MultipatchSplineBuilderXY = MultipatchSplineBuilder2D<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplines1OnPatch,
        BSplines2OnPatch,
        Grid1OnPatch,
        Grid2OnPatch,
        BC,
        BC,
        BC,
        BC,
        BC,
        Connectivity,
        ddc::SplineSolver::LAPACK,
        DConstFieldOnPatch,
        Patch1,
        Patch2>;

class MultipatchSplineBuilder2DTest : public ::testing::Test
{
public:
    MultipatchSplineBuilder2DTest() {}

    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        Coord<X<1>> const x1_min(0.0);
        Coord<X<1>> const x1_max(2 * M_PI);
        IdxStep<GridX<1>> x1_ncells(10);

        Coord<Y<1>> const y1_min(0.0);
        Coord<Y<1>> const y1_max(1.0);
        IdxStep<GridY<1>> y1_ncells(12);

        ddc::init_discrete_space<BSplinesX<1>>(
                build_uniform_break_points(x1_min, x1_max, x1_ncells));
        ddc::init_discrete_space<BSplinesY<1>>(
                build_uniform_break_points(y1_min, y1_max, y1_ncells));

        // Patch 2
        Coord<X<2>> const x2_min(2 * M_PI);
        Coord<X<2>> const x2_max(4 * M_PI);
        IdxStep<GridX<2>> x2_ncells(6);

        Coord<Y<2>> const y2_min(0.0);
        Coord<Y<2>> const y2_max(1.0);
        IdxStep<GridY<2>> y2_ncells(8);

        ddc::init_discrete_space<BSplinesX<2>>(
                build_uniform_break_points(x2_min, x2_max, x2_ncells));
        ddc::init_discrete_space<BSplinesY<2>>(
                build_uniform_break_points(y2_min, y2_max, y2_ncells));
    }

    void initialise_2D_functions(DFieldOnPatch<Patch1> function_1, DFieldOnPatch<Patch2> function_2)
    {
        ddc::parallel_for_each(
                get_idx_range(function_1),
                KOKKOS_LAMBDA(Idx<GridX<1>, GridY<1>> idx) {
                    double const x = ddc::coordinate(Idx<GridX<1>>(idx));
                    double const y = ddc::coordinate(Idx<GridY<1>>(idx));
                    function_1(idx) = Kokkos::sin(x) + Kokkos::cos(y);
                });
        ddc::parallel_for_each(
                get_idx_range(function_2),
                KOKKOS_LAMBDA(Idx<GridX<2>, GridY<2>> idx) {
                    double const x = ddc::coordinate(Idx<GridX<2>>(idx));
                    double const y = ddc::coordinate(Idx<GridY<2>>(idx));
                    function_2(idx) = Kokkos::sin(2 * x * y);
                });
    }

    void initialise_2D_derivatives_2(
            DField<IdxRange<ddc::Deriv<X<2>>, GridY<2>>> derivs_xmin2,
            DField<IdxRange<ddc::Deriv<X<2>>, GridY<2>>> derivs_xmax2,
            DField<IdxRange<GridX<2>, ddc::Deriv<Y<2>>>> derivs_ymin2,
            DField<IdxRange<GridX<2>, ddc::Deriv<Y<2>>>> derivs_ymax2)
    {
        IdxRange<GridX<2>> idx_range_x2 = get_idx_range<GridX<2>>(derivs_ymin2);
        IdxRange<GridY<2>> idx_range_y2 = get_idx_range<GridY<2>>(derivs_xmin2);
        Idx<ddc::Deriv<X<2>>> first_deriv_x2(1);
        Idx<ddc::Deriv<Y<2>>> first_deriv_y2(1);
        double const xmin2 = ddc::coordinate(idx_range_x2.front());
        double const xmax2 = ddc::coordinate(idx_range_x2.back());
        double const ymin2 = ddc::coordinate(idx_range_y2.front());
        double const ymax2 = ddc::coordinate(idx_range_y2.back());
        ddc::parallel_for_each(
                idx_range_y2,
                KOKKOS_LAMBDA(Idx<GridY<2>> idx) {
                    double const y = ddc::coordinate(Idx<GridY<2>>(idx));
                    derivs_xmin2(first_deriv_x2, idx) = 2 * y * Kokkos::cos(2 * xmin2 * y);
                    derivs_xmax2(first_deriv_x2, idx) = 2 * y * Kokkos::cos(2 * xmax2 * y);
                });
        ddc::parallel_for_each(
                idx_range_x2,
                KOKKOS_LAMBDA(Idx<GridX<2>> idx) {
                    double const x = ddc::coordinate(Idx<GridX<2>>(idx));
                    derivs_ymin2(first_deriv_y2, idx) = 2 * x * Kokkos::cos(2 * x * ymin2);
                    derivs_ymax2(first_deriv_y2, idx) = 2 * x * Kokkos::cos(2 * x * ymax2);
                });
    }

    // Test method -------------------------------------------------------------------------------
    void check_if_equal_to_expected(
            MultipatchField<SplineCoeffOnPatch_2D, Patch1, Patch2> const& function_coef,
            SplineCoeffOnPatch_2D<Patch1> function_1_coef_expected,
            SplineCoeffOnPatch_2D<Patch2> function_2_coef_expected)
    {
        // --- get the spline representation on CPU
        auto function_coef_1_host
                = ddc::create_mirror_and_copy(function_coef.template get<Patch1>());
        auto function_coef_2_host
                = ddc::create_mirror_and_copy(function_coef.template get<Patch2>());

        // --- get expected spline representations on CPU
        auto function_coef_1_expected_host = ddc::create_mirror_and_copy(function_1_coef_expected);
        auto function_coef_2_expected_host = ddc::create_mirror_and_copy(function_2_coef_expected);

        // --- check error
        double max_abs_error = 0;
        ddc::for_each(
                get_idx_range(function_coef_1_expected_host),
                [&](Idx<Patch1::BSplines1, Patch1::BSplines2> const idx) {
                    double err
                            = abs(function_coef_1_expected_host(idx) - function_coef_1_host(idx));
                    max_abs_error = std::max(max_abs_error, err);
                });
        ddc::for_each(
                get_idx_range(function_coef_2_expected_host),
                [&](Idx<Patch2::BSplines1, Patch2::BSplines2> const idx) {
                    double err
                            = abs(function_coef_2_expected_host(idx) - function_coef_2_host(idx));
                    max_abs_error = std::max(max_abs_error, err);
                });
        std::cout << "Max absolute error: " << max_abs_error << std::endl;
        EXPECT_LE(max_abs_error, 1e-15);
    }
};

} // namespace


/*
    The global domain is split into two 2D patches on (X1, Y1) and (X2, Y2).
 */
TEST_F(MultipatchSplineBuilder2DTest, TwoPatches2D)
{
    ddc::init_discrete_space<GridX<1>>(
            SplineInterpPointsX<1, ddc::BoundCond::GREVILLE>::get_sampling<GridX<1>>());
    ddc::init_discrete_space<GridY<1>>(
            SplineInterpPointsY<1, ddc::BoundCond::GREVILLE>::get_sampling<GridY<1>>());

    ddc::init_discrete_space<GridX<2>>(
            SplineInterpPointsX<2, ddc::BoundCond::GREVILLE>::get_sampling<GridX<2>>());
    ddc::init_discrete_space<GridY<2>>(
            SplineInterpPointsY<2, ddc::BoundCond::GREVILLE>::get_sampling<GridY<2>>());

    IdxRange<GridX<1>> const idx_range_x1
            = SplineInterpPointsX<1, ddc::BoundCond::GREVILLE>::get_domain<GridX<1>>();
    IdxRange<GridY<1>> const idx_range_y1
            = SplineInterpPointsY<1, ddc::BoundCond::GREVILLE>::get_domain<GridY<1>>();

    IdxRange<GridX<2>> const idx_range_x2
            = SplineInterpPointsX<2, ddc::BoundCond::GREVILLE>::get_domain<GridX<2>>();
    IdxRange<GridY<2>> const idx_range_y2
            = SplineInterpPointsY<2, ddc::BoundCond::GREVILLE>::get_domain<GridY<2>>();

    IdxRange<GridX<1>, GridY<1>> const idx_range_xy1(idx_range_x1, idx_range_y1);
    IdxRange<GridX<2>, GridY<2>> const idx_range_xy2(idx_range_x2, idx_range_y2);

    // List of spline builders
    SplineXYBuilder<1, ddc::BoundCond::GREVILLE> builder1(idx_range_xy1);
    SplineXYBuilder<2, ddc::BoundCond::GREVILLE> builder2(idx_range_xy2);

    // Collection of builders for each patch
    MultipatchSplineBuilderXY<ddc::BoundCond::GREVILLE> builder(builder1, builder2);


    // Function tests
    // --- patch 1
    DFieldMemOnPatch<Patch1> function_1_alloc(idx_range_xy1);
    DFieldOnPatch<Patch1> function_1 = get_field(function_1_alloc);

    // --- patch 2
    DFieldMemOnPatch<Patch2> function_2_alloc(idx_range_xy2);
    DFieldOnPatch<Patch2> function_2 = get_field(function_2_alloc);

    initialise_2D_functions(function_1, function_2);

    // --- collection of values of the function on each patch
    MultipatchField<DConstFieldOnPatch, Patch1, Patch2>
            function_values(get_const_field(function_1), get_const_field(function_2));



    // Spline representations
    // --- patch 1
    DFieldMem<typename Patch1::IdxRangeBS12> function_1_coef_alloc(
            builder1.batched_spline_domain());
    SplineCoeffOnPatch_2D<Patch1> function_1_coef = get_field(function_1_coef_alloc);

    DFieldMem<typename Patch1::IdxRangeBS12> function_1_coef_expected_alloc(
            builder1.batched_spline_domain());
    SplineCoeffOnPatch_2D<Patch1> function_1_coef_expected
            = get_field(function_1_coef_expected_alloc);

    // --- patch 2
    DFieldMem<typename Patch2::IdxRangeBS12> function_2_coef_alloc(
            builder2.batched_spline_domain());
    SplineCoeffOnPatch_2D<Patch2> function_2_coef = get_field(function_2_coef_alloc);

    DFieldMem<typename Patch2::IdxRangeBS12> function_2_coef_expected_alloc(
            builder2.batched_spline_domain());
    SplineCoeffOnPatch_2D<Patch2> function_2_coef_expected
            = get_field(function_2_coef_expected_alloc);

    // --- collection
    MultipatchField<SplineCoeffOnPatch_2D, Patch1, Patch2>
            function_coef(function_1_coef, function_2_coef);


    // Build the spline representations on each patch
    builder(function_coef, function_values);



    // Check the results are the same than with one by one spline building
    builder1(function_1_coef_expected, get_const_field(function_1));
    builder2(function_2_coef_expected, get_const_field(function_2));

    check_if_equal_to_expected(function_coef, function_1_coef_expected, function_2_coef_expected);
}

TEST_F(MultipatchSplineBuilder2DTest, TwoPatches2DHermite)
{
    ddc::init_discrete_space<GridX<1>>(
            SplineInterpPointsX<1, ddc::BoundCond::HERMITE>::get_sampling<GridX<1>>());
    ddc::init_discrete_space<GridY<1>>(
            SplineInterpPointsY<1, ddc::BoundCond::HERMITE>::get_sampling<GridY<1>>());

    ddc::init_discrete_space<GridX<2>>(
            SplineInterpPointsX<2, ddc::BoundCond::HERMITE>::get_sampling<GridX<2>>());
    ddc::init_discrete_space<GridY<2>>(
            SplineInterpPointsY<2, ddc::BoundCond::HERMITE>::get_sampling<GridY<2>>());

    IdxRange<GridX<1>> const idx_range_x1
            = SplineInterpPointsX<1, ddc::BoundCond::HERMITE>::get_domain<GridX<1>>();
    IdxRange<GridY<1>> const idx_range_y1
            = SplineInterpPointsY<1, ddc::BoundCond::HERMITE>::get_domain<GridY<1>>();

    IdxRange<GridX<2>> const idx_range_x2
            = SplineInterpPointsX<2, ddc::BoundCond::HERMITE>::get_domain<GridX<2>>();
    IdxRange<GridY<2>> const idx_range_y2
            = SplineInterpPointsY<2, ddc::BoundCond::HERMITE>::get_domain<GridY<2>>();

    IdxRange<GridX<1>, GridY<1>> const idx_range_xy1(idx_range_x1, idx_range_y1);
    IdxRange<GridX<2>, GridY<2>> const idx_range_xy2(idx_range_x2, idx_range_y2);

    // List of spline builders
    SplineXYBuilder<1, ddc::BoundCond::HERMITE> builder1(idx_range_xy1);
    SplineXYBuilder<2, ddc::BoundCond::HERMITE> builder2(idx_range_xy2);

    // Collection of builders for each patch
    MultipatchSplineBuilderXY<ddc::BoundCond::HERMITE> builder(builder1, builder2);


    // Function tests
    // --- patch 1
    DFieldMemOnPatch<Patch1> function_1_alloc(idx_range_xy1);
    DFieldOnPatch<Patch1> function_1 = get_field(function_1_alloc);

    // --- patch 2
    DFieldMemOnPatch<Patch2> function_2_alloc(idx_range_xy2);
    DFieldOnPatch<Patch2> function_2 = get_field(function_2_alloc);

    initialise_2D_functions(function_1, function_2);

    // --- collection of values of the function on each patch
    MultipatchField<DConstFieldOnPatch, Patch1, Patch2>
            function_values(get_const_field(function_1), get_const_field(function_2));

    // Derivatives Patch1
    Idx<ddc::Deriv<X<1>>> first_deriv_x1(1);
    IdxStep<ddc::Deriv<X<1>>> n_deriv_x1(1);
    IdxRange<ddc::Deriv<X<1>>> deriv_x_idx_range1(first_deriv_x1, n_deriv_x1);

    Idx<ddc::Deriv<Y<1>>> first_deriv_y1(1);
    IdxStep<ddc::Deriv<Y<1>>> n_deriv_y1(1);
    IdxRange<ddc::Deriv<Y<1>>> deriv_y_idx_range1(first_deriv_y1, n_deriv_y1);

    IdxRange<ddc::Deriv<X<1>>, GridY<1>> derivs_x1_idx_range(deriv_x_idx_range1, idx_range_y1);
    DFieldMem<IdxRange<ddc::Deriv<X<1>>, GridY<1>>> derivs_xmin1_alloc(derivs_x1_idx_range);
    DField<IdxRange<ddc::Deriv<X<1>>, GridY<1>>> derivs_xmin1 = get_field(derivs_xmin1_alloc);
    DFieldMem<IdxRange<ddc::Deriv<X<1>>, GridY<1>>> derivs_xmax1_alloc(derivs_x1_idx_range);
    DField<IdxRange<ddc::Deriv<X<1>>, GridY<1>>> derivs_xmax1 = get_field(derivs_xmax1_alloc);

    double const xmin1 = ddc::coordinate(idx_range_x1.front());
    double const xmax1 = ddc::coordinate(idx_range_x1.back());
    ddc::parallel_fill(derivs_xmin1, Kokkos::cos(xmin1));
    ddc::parallel_fill(derivs_xmax1, Kokkos::cos(xmax1));

    IdxRange<GridX<1>, ddc::Deriv<Y<1>>> derivs_y1_idx_range(deriv_y_idx_range1, idx_range_x1);
    DFieldMem<IdxRange<GridX<1>, ddc::Deriv<Y<1>>>> derivs_ymin1_alloc(derivs_y1_idx_range);
    DField<IdxRange<GridX<1>, ddc::Deriv<Y<1>>>> derivs_ymin1 = get_field(derivs_ymin1_alloc);
    DFieldMem<IdxRange<GridX<1>, ddc::Deriv<Y<1>>>> derivs_ymax1_alloc(derivs_y1_idx_range);
    DField<IdxRange<GridX<1>, ddc::Deriv<Y<1>>>> derivs_ymax1 = get_field(derivs_ymax1_alloc);

    double const ymin1 = ddc::coordinate(idx_range_y1.front());
    double const ymax1 = ddc::coordinate(idx_range_y1.back());
    ddc::parallel_fill(derivs_ymin1, -Kokkos::sin(ymin1));
    ddc::parallel_fill(derivs_ymax1, -Kokkos::sin(ymax1));

    // Derivatives Patch 2
    Idx<ddc::Deriv<X<2>>> first_deriv_x2(1);
    IdxStep<ddc::Deriv<X<2>>> n_deriv_x2(1);
    IdxRange<ddc::Deriv<X<2>>> deriv_x_idx_range2(first_deriv_x2, n_deriv_x2);

    Idx<ddc::Deriv<Y<2>>> first_deriv_y2(1);
    IdxStep<ddc::Deriv<Y<2>>> n_deriv_y2(1);
    IdxRange<ddc::Deriv<Y<2>>> deriv_y_idx_range2(first_deriv_y2, n_deriv_y2);

    IdxRange<ddc::Deriv<X<2>>, GridY<2>> derivs_x2_idx_range(deriv_x_idx_range2, idx_range_y2);
    DFieldMem<IdxRange<ddc::Deriv<X<2>>, GridY<2>>> derivs_xmin2_alloc(derivs_x2_idx_range);
    DField<IdxRange<ddc::Deriv<X<2>>, GridY<2>>> derivs_xmin2 = get_field(derivs_xmin2_alloc);
    DFieldMem<IdxRange<ddc::Deriv<X<2>>, GridY<2>>> derivs_xmax2_alloc(derivs_x2_idx_range);
    DField<IdxRange<ddc::Deriv<X<2>>, GridY<2>>> derivs_xmax2 = get_field(derivs_xmax2_alloc);

    IdxRange<GridX<2>, ddc::Deriv<Y<2>>> derivs_y2_idx_range(deriv_y_idx_range2, idx_range_x2);
    DFieldMem<IdxRange<GridX<2>, ddc::Deriv<Y<2>>>> derivs_ymin2_alloc(derivs_y2_idx_range);
    DField<IdxRange<GridX<2>, ddc::Deriv<Y<2>>>> derivs_ymin2 = get_field(derivs_ymin2_alloc);
    DFieldMem<IdxRange<GridX<2>, ddc::Deriv<Y<2>>>> derivs_ymax2_alloc(derivs_y2_idx_range);
    DField<IdxRange<GridX<2>, ddc::Deriv<Y<2>>>> derivs_ymax2 = get_field(derivs_ymax2_alloc);

    initialise_2D_derivatives_2(derivs_xmin2, derivs_xmax2, derivs_ymin2, derivs_ymax2);


    // --- collection of values of the function on each patch
    MultipatchField<ConstDeriv1_OnPatch_2D, Patch1, Patch2>
            derivs_xmin(get_const_field(derivs_xmin1), get_const_field(derivs_xmin2));
    MultipatchField<ConstDeriv1_OnPatch_2D, Patch1, Patch2>
            derivs_xmax(get_const_field(derivs_xmax1), get_const_field(derivs_xmax2));
    MultipatchField<ConstDeriv2_OnPatch_2D, Patch1, Patch2>
            derivs_ymin(get_const_field(derivs_ymin1), get_const_field(derivs_ymin2));
    MultipatchField<ConstDeriv2_OnPatch_2D, Patch1, Patch2>
            derivs_ymax(get_const_field(derivs_ymax1), get_const_field(derivs_ymax2));

    IdxRange<ddc::Deriv<X<1>>, ddc::Deriv<Y<1>>>
            cross_derivs_idx_range1(deriv_x_idx_range1, deriv_y_idx_range1);
    IdxRange<ddc::Deriv<X<2>>, ddc::Deriv<Y<2>>>
            cross_derivs_idx_range2(deriv_x_idx_range2, deriv_y_idx_range2);

    DFieldMem<IdxRange<ddc::Deriv<X<1>>, ddc::Deriv<Y<1>>>> mixed_derivs_min1_min2_alloc1(
            cross_derivs_idx_range1);
    DFieldMem<IdxRange<ddc::Deriv<X<1>>, ddc::Deriv<Y<1>>>> mixed_derivs_max1_min2_alloc1(
            cross_derivs_idx_range1);
    DFieldMem<IdxRange<ddc::Deriv<X<1>>, ddc::Deriv<Y<1>>>> mixed_derivs_min1_max2_alloc1(
            cross_derivs_idx_range1);
    DFieldMem<IdxRange<ddc::Deriv<X<1>>, ddc::Deriv<Y<1>>>> mixed_derivs_max1_max2_alloc1(
            cross_derivs_idx_range1);

    DFieldMem<IdxRange<ddc::Deriv<X<2>>, ddc::Deriv<Y<2>>>> mixed_derivs_min1_min2_alloc2(
            cross_derivs_idx_range2);
    DFieldMem<IdxRange<ddc::Deriv<X<2>>, ddc::Deriv<Y<2>>>> mixed_derivs_max1_min2_alloc2(
            cross_derivs_idx_range2);
    DFieldMem<IdxRange<ddc::Deriv<X<2>>, ddc::Deriv<Y<2>>>> mixed_derivs_min1_max2_alloc2(
            cross_derivs_idx_range2);
    DFieldMem<IdxRange<ddc::Deriv<X<2>>, ddc::Deriv<Y<2>>>> mixed_derivs_max1_max2_alloc2(
            cross_derivs_idx_range2);

    ddc::parallel_fill(get_field(mixed_derivs_min1_min2_alloc1), 0.0);
    ddc::parallel_fill(get_field(mixed_derivs_max1_min2_alloc1), 0.0);
    ddc::parallel_fill(get_field(mixed_derivs_min1_max2_alloc1), 0.0);
    ddc::parallel_fill(get_field(mixed_derivs_max1_max2_alloc1), 0.0);

    double const xmin2 = ddc::coordinate(idx_range_x2.front());
    double const xmax2 = ddc::coordinate(idx_range_x2.back());
    double const ymin2 = ddc::coordinate(idx_range_y2.front());
    double const ymax2 = ddc::coordinate(idx_range_y2.back());
    ddc::parallel_fill(
            get_field(mixed_derivs_min1_min2_alloc2),
            2.
                    * (Kokkos::cos(2. * xmin2 * ymin2)
                       - 2. * xmin2 * ymin2 * Kokkos::sin(2. * xmin2 * ymin2)));
    ddc::parallel_fill(
            get_field(mixed_derivs_max1_min2_alloc2),
            2.
                    * (Kokkos::cos(2. * xmax2 * ymin2)
                       - 2. * xmax2 * ymin2 * Kokkos::sin(2. * xmax2 * ymin2)));
    ddc::parallel_fill(
            get_field(mixed_derivs_min1_max2_alloc2),
            2.
                    * (Kokkos::cos(2. * xmin2 * ymax2)
                       - 2. * xmin2 * ymax2 * Kokkos::sin(2. * xmin2 * ymax2)));
    ddc::parallel_fill(
            get_field(mixed_derivs_max1_max2_alloc2),
            2.
                    * (Kokkos::cos(2. * xmax2 * ymax2)
                       - 2. * xmax2 * ymax2 * Kokkos::sin(2. * xmax2 * ymax2)));

    MultipatchField<ConstDeriv12_OnPatch_2D, Patch1, Patch2> mixed_derivs_min1_min2(
            get_const_field(mixed_derivs_min1_min2_alloc1),
            get_const_field(mixed_derivs_min1_min2_alloc2));
    MultipatchField<ConstDeriv12_OnPatch_2D, Patch1, Patch2> mixed_derivs_max1_min2(
            get_const_field(mixed_derivs_max1_min2_alloc1),
            get_const_field(mixed_derivs_max1_min2_alloc2));
    MultipatchField<ConstDeriv12_OnPatch_2D, Patch1, Patch2> mixed_derivs_min1_max2(
            get_const_field(mixed_derivs_min1_max2_alloc1),
            get_const_field(mixed_derivs_min1_max2_alloc2));
    MultipatchField<ConstDeriv12_OnPatch_2D, Patch1, Patch2> mixed_derivs_max1_max2(
            get_const_field(mixed_derivs_max1_max2_alloc1),
            get_const_field(mixed_derivs_max1_max2_alloc2));

    // Spline representations
    // --- patch 1
    DFieldMem<typename Patch1::IdxRangeBS12> function_1_coef_alloc(
            builder1.batched_spline_domain());
    SplineCoeffOnPatch_2D<Patch1> function_1_coef = get_field(function_1_coef_alloc);

    DFieldMem<typename Patch1::IdxRangeBS12> function_1_coef_expected_alloc(
            builder1.batched_spline_domain());
    SplineCoeffOnPatch_2D<Patch1> function_1_coef_expected
            = get_field(function_1_coef_expected_alloc);

    // --- patch 2
    DFieldMem<typename Patch2::IdxRangeBS12> function_2_coef_alloc(
            builder2.batched_spline_domain());
    SplineCoeffOnPatch_2D<Patch2> function_2_coef = get_field(function_2_coef_alloc);

    DFieldMem<typename Patch2::IdxRangeBS12> function_2_coef_expected_alloc(
            builder2.batched_spline_domain());
    SplineCoeffOnPatch_2D<Patch2> function_2_coef_expected
            = get_field(function_2_coef_expected_alloc);

    // --- collection
    MultipatchField<SplineCoeffOnPatch_2D, Patch1, Patch2>
            function_coef(function_1_coef, function_2_coef);


    // Build the spline representations on each patch
    builder(function_coef,
            function_values,
            derivs_xmin,
            derivs_xmax,
            derivs_ymin,
            derivs_ymax,
            mixed_derivs_min1_min2,
            mixed_derivs_max1_min2,
            mixed_derivs_min1_max2,
            mixed_derivs_max1_max2);



    // Check the results are the same than with one by one spline building
    builder1(
            function_1_coef_expected,
            get_const_field(function_1),
            std::optional(get_const_field(derivs_xmin1)),
            std::optional(get_const_field(derivs_xmax1)),
            std::optional(get_const_field(derivs_ymin1)),
            std::optional(get_const_field(derivs_ymax1)),
            std::optional(get_const_field(mixed_derivs_min1_min2_alloc1)),
            std::optional(get_const_field(mixed_derivs_max1_min2_alloc1)),
            std::optional(get_const_field(mixed_derivs_min1_max2_alloc1)),
            std::optional(get_const_field(mixed_derivs_max1_max2_alloc1)));
    builder2(
            function_2_coef_expected,
            get_const_field(function_2),
            std::optional(get_const_field(derivs_xmin2)),
            std::optional(get_const_field(derivs_xmax2)),
            std::optional(get_const_field(derivs_ymin2)),
            std::optional(get_const_field(derivs_ymax2)),
            std::optional(get_const_field(mixed_derivs_min1_min2_alloc2)),
            std::optional(get_const_field(mixed_derivs_max1_min2_alloc2)),
            std::optional(get_const_field(mixed_derivs_min1_max2_alloc2)),
            std::optional(get_const_field(mixed_derivs_max1_max2_alloc2)));

    check_if_equal_to_expected(function_coef, function_1_coef_expected, function_2_coef_expected);
}
