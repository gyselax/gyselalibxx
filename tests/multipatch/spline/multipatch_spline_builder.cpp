// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_non_periodic_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "mesh_builder.hpp"
#include "multipatch_field.hpp"
#include "multipatch_spline_builder.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"


namespace {
using namespace non_periodic_non_uniform_2d_2patches;
ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::GREVILLE;
ddc::BoundCond constexpr SplineYBoundary = ddc::BoundCond::GREVILLE;

ddc::BoundCond constexpr SplineX1Boundary = SplineXBoundary;
ddc::BoundCond constexpr SplineY1Boundary = SplineYBoundary;

ddc::BoundCond constexpr SplineX2Boundary = SplineXBoundary;
ddc::BoundCond constexpr SplineY2Boundary = SplineYBoundary;

using SplineInterpPointsX1
        = ddc::GrevilleInterpolationPoints<BSplinesX<1>, SplineX1Boundary, SplineX1Boundary>;
using SplineInterpPointsY1
        = ddc::GrevilleInterpolationPoints<BSplinesY<1>, SplineY1Boundary, SplineY1Boundary>;
using SplineInterpPointsX2

        = ddc::GrevilleInterpolationPoints<BSplinesX<2>, SplineX2Boundary, SplineX2Boundary>;
using SplineInterpPointsY2
        = ddc::GrevilleInterpolationPoints<BSplinesY<2>, SplineY2Boundary, SplineY2Boundary>;

// Operators
using SplineX1Builder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX<1>,
        GridX<1>,
        SplineX1Boundary,
        SplineX1Boundary,
        ddc::SplineSolver::LAPACK>;

using SplineX2Builder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX<2>,
        GridX<2>,
        SplineX2Boundary,
        SplineX2Boundary,
        ddc::SplineSolver::LAPACK>;



using MultipatchSplineBuilderX = MultipatchSplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplines1OnPatch,
        Grid1OnPatch,
        SplineXBoundary,
        SplineXBoundary,
        SplineXBoundary,
        Connectivity,
        ddc::SplineSolver::LAPACK,
        DConstFieldOnPatch,
        Patch1,
        Patch2>;


class MultipatchSplineBuilderTest : public ::testing::Test
{
private:
    static constexpr IdxStep<GridX<1>> x1_size = IdxStep<GridX<1>>(10);
    static constexpr IdxStep<GridY<1>> y1_size = IdxStep<GridY<1>>(10);

    static constexpr IdxStep<GridX<2>> x2_size = IdxStep<GridX<2>>(10);
    static constexpr IdxStep<GridY<2>> y2_size = IdxStep<GridY<2>>(10);

protected:
    IdxRange<GridX<1>> const idx_range_x1;
    IdxRange<GridY<1>> const idx_range_y1;

    IdxRange<GridX<2>> const idx_range_x2;
    IdxRange<GridY<2>> const idx_range_y2;

public:
    MultipatchSplineBuilderTest()
        : idx_range_x1(SplineInterpPointsX1::get_domain<GridX<1>>())
        , idx_range_y1(SplineInterpPointsY1::get_domain<GridY<1>>())
        , idx_range_x2(SplineInterpPointsX2::get_domain<GridX<2>>())
        , idx_range_y2(SplineInterpPointsY2::get_domain<GridY<2>>()) {};

    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        Coord<X<1>> const x1_min(0.0);
        Coord<X<1>> const x1_max(2 * M_PI);

        Coord<Y<1>> const y1_min(0.0);
        Coord<Y<1>> const y1_max(1.0);

        ddc::init_discrete_space<BSplinesX<1>>(build_uniform_break_points(x1_min, x1_max, x1_size));
        ddc::init_discrete_space<BSplinesY<1>>(build_uniform_break_points(y1_min, y1_max, y1_size));

        ddc::init_discrete_space<GridX<1>>(SplineInterpPointsX1::get_sampling<GridX<1>>());
        ddc::init_discrete_space<GridY<1>>(SplineInterpPointsY1::get_sampling<GridY<1>>());

        // Patch 2
        Coord<X<2>> const x2_min(2 * M_PI);
        Coord<X<2>> const x2_max(4 * M_PI);

        Coord<Y<2>> const y2_min(0.0);
        Coord<Y<2>> const y2_max(1.0);

        ddc::init_discrete_space<BSplinesX<2>>(build_uniform_break_points(x2_min, x2_max, x2_size));
        ddc::init_discrete_space<BSplinesY<2>>(build_uniform_break_points(y2_min, y2_max, y2_size));

        ddc::init_discrete_space<GridX<2>>(SplineInterpPointsX2::get_sampling<GridX<2>>());
        ddc::init_discrete_space<GridY<2>>(SplineInterpPointsY2::get_sampling<GridY<2>>());
    }


    // Initialisation methods --------------------------------------------------------------------
    void initialise_1D_functions(
            Field<double, IdxRange<GridX<1>>> function_1,
            Field<double, IdxRange<GridX<2>>> function_2)
    {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_x1,
                KOKKOS_LAMBDA(Idx<GridX<1>> const idx) {
                    double const x = ddc::coordinate(idx);
                    function_1(idx) = Kokkos::sin(x);
                });

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_x2,
                KOKKOS_LAMBDA(Idx<GridX<2>> const idx) {
                    double const x = ddc::coordinate(idx);
                    function_2(idx) = Kokkos::sin(2 * x);
                });
    }


    void initialise_2D_functions(
            Field<double, IdxRange<GridX<1>, GridY<1>>> function_1,
            Field<double, IdxRange<GridX<2>, GridY<2>>> function_2)
    {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(function_1),
                KOKKOS_LAMBDA(Idx<GridX<1>, GridY<1>> const idx) {
                    double const x = ddc::coordinate(Idx<GridX<1>>(idx));
                    double const y = ddc::coordinate(Idx<GridY<1>>(idx));
                    function_1(idx) = Kokkos::sin(x) + Kokkos::cos(y);
                });

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(function_2),
                KOKKOS_LAMBDA(Idx<GridX<2>, GridY<2>> const idx) {
                    double const x = ddc::coordinate(Idx<GridX<2>>(idx));
                    double const y = ddc::coordinate(Idx<GridY<2>>(idx));
                    function_2(idx) = Kokkos::sin(2 * x) + Kokkos::cos(y);
                });
    }


    // Test method -------------------------------------------------------------------------------
    template <class MultipatchCoeff, class SplineField1, class SplineField2>
    void check_if_equal_to_expected(
            MultipatchCoeff const& function_coef,
            SplineField1 const& function_1_coef_expected,
            SplineField2 const& function_2_coef_expected)
    {
        using Idx1 = typename SplineField1::discrete_domain_type::discrete_element_type;
        using Idx2 = typename SplineField2::discrete_domain_type::discrete_element_type;

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
        ddc::for_each(get_idx_range(function_coef_1_expected_host), [&](Idx1 const idx) {
            double err = abs(function_coef_1_expected_host(idx) - function_coef_1_host(idx));
            max_abs_error = std::max(max_abs_error, err);
        });
        ddc::for_each(get_idx_range(function_coef_2_expected_host), [&](Idx2 const idx) {
            double err = abs(function_coef_2_expected_host(idx) - function_coef_2_host(idx));
            max_abs_error = std::max(max_abs_error, err);
        });
        EXPECT_LE(max_abs_error, 1e-15);
        std::cout << "Max absolute error: " << max_abs_error << std::endl;
    }
};


} // namespace


/*
    The global domain is split into two 1D patches on X1 and X2.
 */
TEST_F(MultipatchSplineBuilderTest, TwoPatches1D)
{
    // List of spline builders
    SplineX1Builder builder_x1(idx_range_x1);
    SplineX2Builder builder_x2(idx_range_x2);

    // Collection of builders for each patch
    MultipatchSplineBuilderX builder(builder_x1, builder_x2);


    // Function tests
    // --- patch 1
    DFieldMem<IdxRange<GridX<1>>> function_1_alloc(idx_range_x1);
    DField<IdxRange<GridX<1>>> function_1 = get_field(function_1_alloc);

    // --- patch 2
    DFieldMem<IdxRange<GridX<2>>> function_2_alloc(idx_range_x2);
    DField<IdxRange<GridX<2>>> function_2 = get_field(function_2_alloc);

    initialise_1D_functions(function_1, function_2);

    // --- collection of values of the function on each patch
    MultipatchField<DConstField1OnPatch, Patch1, Patch2>
            function_values(get_const_field(function_1), get_const_field(function_2));



    // Spline representations
    // --- patch 1
    DFieldMem<IdxRange<BSplinesX<1>>> function_1_coef_alloc(get_spline_idx_range(builder_x1));
    DField<IdxRange<BSplinesX<1>>> function_1_coef = get_field(function_1_coef_alloc);

    DFieldMem<IdxRange<BSplinesX<1>>> function_1_coef_expected_alloc(
            get_spline_idx_range(builder_x1));
    DField<IdxRange<BSplinesX<1>>> function_1_coef_expected
            = get_field(function_1_coef_expected_alloc);

    // --- patch 2
    DFieldMem<IdxRange<BSplinesX<2>>> function_2_coef_alloc(get_spline_idx_range(builder_x2));
    DField<IdxRange<BSplinesX<2>>> function_2_coef = get_field(function_2_coef_alloc);

    DFieldMem<IdxRange<BSplinesX<2>>> function_2_coef_expected_alloc(
            get_spline_idx_range(builder_x2));
    DField<IdxRange<BSplinesX<2>>> function_2_coef_expected
            = get_field(function_2_coef_expected_alloc);

    // --- collection
    MultipatchField<SplineCoeff1OnPatch_1D, Patch1, Patch2>
            function_coef(function_1_coef, function_2_coef);


    // Build the spline representations on each patch
    builder(function_coef, function_values);



    // Check the results are the same than with one by one spline building
    builder_x1(function_1_coef_expected, get_const_field(function_1));
    builder_x2(function_2_coef_expected, get_const_field(function_2));
    check_if_equal_to_expected(function_coef, function_1_coef_expected, function_2_coef_expected);
}



/*
    The global domain is split into two 2D patches on (X1, Y1) and (X2, Y2).
 */
TEST_F(MultipatchSplineBuilderTest, TwoPatches2D)
{
    // List of index ranges
    IdxRange<GridX<1>, GridY<1>> idx_range_xy1(idx_range_x1, idx_range_y1);
    IdxRange<GridX<2>, GridY<2>> idx_range_xy2(idx_range_x2, idx_range_y2);


    // List of spline builders
    SplineX1Builder builder_x1(idx_range_xy1);
    SplineX2Builder builder_x2(idx_range_xy2);

    // Collection of builders for each patch
    MultipatchSplineBuilderX builder(builder_x1, builder_x2);


    // Function tests
    // --- patch 1
    DFieldMem<IdxRange<GridX<1>, GridY<1>>> function_1_alloc(idx_range_xy1);
    DField<IdxRange<GridX<1>, GridY<1>>> function_1 = get_field(function_1_alloc);

    // --- patch 2
    DFieldMem<IdxRange<GridX<2>, GridY<2>>> function_2_alloc(idx_range_xy2);
    DField<IdxRange<GridX<2>, GridY<2>>> function_2 = get_field(function_2_alloc);

    initialise_2D_functions(function_1, function_2);

    // --- collection of values of the function on each patch
    MultipatchField<DConstFieldOnPatch, Patch1, Patch2>
            function_values(get_const_field(function_1), get_const_field(function_2));



    // Spline representations
    // --- patch 1
    DFieldMem<IdxRange<BSplinesX<1>, GridY<1>>> function_1_coef_alloc(
            builder_x1.batched_spline_domain());
    DField<IdxRange<BSplinesX<1>, GridY<1>>> function_1_coef = get_field(function_1_coef_alloc);

    DFieldMem<IdxRange<BSplinesX<1>, GridY<1>>> function_1_coef_expected_alloc(
            builder_x1.batched_spline_domain());
    DField<IdxRange<BSplinesX<1>, GridY<1>>> function_1_coef_expected
            = get_field(function_1_coef_expected_alloc);

    // --- patch 2
    DFieldMem<IdxRange<BSplinesX<2>, GridY<2>>> function_2_coef_alloc(
            builder_x2.batched_spline_domain());
    DField<IdxRange<BSplinesX<2>, GridY<2>>> function_2_coef = get_field(function_2_coef_alloc);

    DFieldMem<IdxRange<BSplinesX<2>, GridY<2>>> function_2_coef_expected_alloc(
            builder_x2.batched_spline_domain());
    DField<IdxRange<BSplinesX<2>, GridY<2>>> function_2_coef_expected
            = get_field(function_2_coef_expected_alloc);

    // --- collection
    MultipatchField<SplineCoeff1OnPatch_2D, Patch1, Patch2>
            function_coef(function_1_coef, function_2_coef);


    // Build the spline representations on each patch
    builder(function_coef, function_values);



    // Check the results are the same than with one by one spline building
    builder_x1(function_1_coef_expected, get_const_field(function_1));
    builder_x2(function_2_coef_expected, get_const_field(function_2));
    check_if_equal_to_expected(function_coef, function_1_coef_expected, function_2_coef_expected);
}
