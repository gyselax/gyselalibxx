// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "ddc_helper.hpp"
#include "directional_tag.hpp"
#include "multipatch_spline_builder.hpp"
#include "vector_field.hpp"
#include "vector_field_span.hpp"


namespace {
// Continuous dimension of patch 1
struct X1
{
    static bool constexpr PERIODIC = false;
};

struct Y1
{
    static bool constexpr PERIODIC = false;
};

// Continuous dimension of patch 2
struct X2
{
    static bool constexpr PERIODIC = false;
};

struct Y2
{
    static bool constexpr PERIODIC = false;
};



// Splines
struct BSplinesX1 : ddc::UniformBSplines<X1, 3>
{
};
struct BSplinesY1 : ddc::UniformBSplines<Y1, 3>
{
};
struct BSplinesX2 : ddc::UniformBSplines<X2, 3>
{
};
struct BSplinesY2 : ddc::UniformBSplines<Y2, 3>
{
};



ddc::BoundCond constexpr SplineX1Boundary = ddc::BoundCond::GREVILLE;
ddc::BoundCond constexpr SplineY1Boundary = ddc::BoundCond::GREVILLE;

ddc::BoundCond constexpr SplineX2Boundary = ddc::BoundCond::GREVILLE;
ddc::BoundCond constexpr SplineY2Boundary = ddc::BoundCond::GREVILLE;



using SplineInterpPointsX1
        = ddc::GrevilleInterpolationPoints<BSplinesX1, SplineX1Boundary, SplineX1Boundary>;
using SplineInterpPointsY1
        = ddc::GrevilleInterpolationPoints<BSplinesY1, SplineY1Boundary, SplineY1Boundary>;


using SplineInterpPointsX2
        = ddc::GrevilleInterpolationPoints<BSplinesX2, SplineX2Boundary, SplineX2Boundary>;
using SplineInterpPointsY2
        = ddc::GrevilleInterpolationPoints<BSplinesY2, SplineY2Boundary, SplineY2Boundary>;


// Discrete dimension of patch 1
struct GridX1 : SplineInterpPointsX1::interpolation_discrete_dimension_type
{
};
struct GridY1 : SplineInterpPointsY1::interpolation_discrete_dimension_type
{
};
// Discrete dimension of patch 2
struct GridX2 : SplineInterpPointsX2::interpolation_discrete_dimension_type
{
};
struct GridY2 : SplineInterpPointsY2::interpolation_discrete_dimension_type
{
};



/*
    TODO: use Patch class once merged. 
*/
using CoordX1 = Coord<X1>;
using CoordY1 = Coord<Y1>;
using CoordX2 = Coord<X2>;
using CoordY2 = Coord<Y2>;


using IdxX1 = Idx<GridX1>;
using IdxY1 = Idx<GridY1>;
using IdxX2 = Idx<GridX2>;
using IdxY2 = Idx<GridY2>;

using IdxXY1 = Idx<GridX1, GridY1>;
using IdxXY2 = Idx<GridX2, GridY2>;


using IdxStepX1 = IdxStep<GridX1>;
using IdxStepY1 = IdxStep<GridY1>;
using IdxStepX2 = IdxStep<GridX2>;
using IdxStepY2 = IdxStep<GridY2>;


using IdxRangeX1 = IdxRange<GridX1>;
using IdxRangeY1 = IdxRange<GridY1>;
using IdxRangeX2 = IdxRange<GridX2>;
using IdxRangeY2 = IdxRange<GridY2>;

using IdxRangeXY1 = IdxRange<GridX1, GridY1>;
using IdxRangeXY2 = IdxRange<GridX2, GridY2>;


// Operators
using SplineX1Builder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX1,
        GridX1,
        SplineX1Boundary,
        SplineX1Boundary,
        ddc::SplineSolver::LAPACK,
        GridX1>;

using SplineX2Builder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX2,
        GridX2,
        SplineX2Boundary,
        SplineX2Boundary,
        ddc::SplineSolver::LAPACK,
        GridX2>;



using SplineX1Builder_2d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX1,
        GridX1,
        SplineX1Boundary,
        SplineX1Boundary,
        ddc::SplineSolver::LAPACK,
        GridX1,
        GridY1>;

using SplineX2Builder_2d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX2,
        GridX2,
        SplineX2Boundary,
        SplineX2Boundary,
        ddc::SplineSolver::LAPACK,
        GridX2,
        GridY2>;



class MultipatchSplineBuilderTest : public ::testing::Test
{
private:
    static constexpr IdxStepX1 x1_size = IdxStepX1(10);
    static constexpr IdxStepY1 y1_size = IdxStepY1(10);

    static constexpr IdxStepX2 x2_size = IdxStepX2(10);
    static constexpr IdxStepY2 y2_size = IdxStepY2(10);

protected:
    IdxRangeX1 const idx_range_x1;
    IdxRangeY1 const idx_range_y1;

    IdxRangeX2 const idx_range_x2;
    IdxRangeY2 const idx_range_y2;

public:
    MultipatchSplineBuilderTest()
        : idx_range_x1(SplineInterpPointsX1::get_domain<GridX1>())
        , idx_range_y1(SplineInterpPointsY1::get_domain<GridY1>())
        , idx_range_x2(SplineInterpPointsX2::get_domain<GridX2>())
        , idx_range_y2(SplineInterpPointsY2::get_domain<GridY2>()) {};

    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        CoordX1 const x1_min(0.0);
        CoordX1 const x1_max(2 * M_PI);

        CoordY1 const y1_min(0.0);
        CoordY1 const y1_max(1.0);

        ddc::init_discrete_space<BSplinesX1>(x1_min, x1_max, x1_size);
        ddc::init_discrete_space<BSplinesY1>(y1_min, y1_max, y1_size);

        ddc::init_discrete_space<GridX1>(SplineInterpPointsX1::get_sampling<GridX1>());
        ddc::init_discrete_space<GridY1>(SplineInterpPointsY1::get_sampling<GridY1>());

        // Patch 2
        CoordX2 const x2_min(2 * M_PI);
        CoordX2 const x2_max(4 * M_PI);

        CoordY2 const y2_min(0.0);
        CoordY2 const y2_max(1.0);

        ddc::init_discrete_space<BSplinesX2>(x2_min, x2_max, x2_size);
        ddc::init_discrete_space<BSplinesY2>(y2_min, y2_max, y2_size);

        ddc::init_discrete_space<GridX2>(SplineInterpPointsX2::get_sampling<GridX2>());
        ddc::init_discrete_space<GridY2>(SplineInterpPointsY2::get_sampling<GridY2>());
    }


    // Initialisation methods --------------------------------------------------------------------
    void initialize_1D_functions(
            Field<double, IdxRangeX1> function_1,
            Field<double, IdxRangeX2> function_2)
    {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_x1,
                KOKKOS_LAMBDA(IdxX1 const idx) {
                    double const x = ddc::coordinate(idx);
                    function_1(idx) = Kokkos::sin(x);
                });

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_x2,
                KOKKOS_LAMBDA(IdxX2 const idx) {
                    double const x = ddc::coordinate(idx);
                    function_2(idx) = Kokkos::sin(2 * x);
                });
    }


    void initialize_2D_functions(
            Field<double, IdxRangeXY1> function_1,
            Field<double, IdxRangeXY2> function_2)
    {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(function_1),
                KOKKOS_LAMBDA(IdxXY1 const idx) {
                    double const x = ddc::coordinate(IdxX1(idx));
                    double const y = ddc::coordinate(IdxY1(idx));
                    function_1(idx) = Kokkos::sin(x) + Kokkos::cos(y);
                });

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(function_2),
                KOKKOS_LAMBDA(IdxXY2 const idx) {
                    double const x = ddc::coordinate(IdxX2(idx));
                    double const y = ddc::coordinate(IdxY2(idx));
                    function_2(idx) = Kokkos::sin(2 * x) + Kokkos::cos(y);
                });
    }


    // Test method -------------------------------------------------------------------------------
    template <class CoefTuple, class SplineField1, class SplineField2>
    void check_if_equal_to_expected(
            CoefTuple const& function_coef,
            SplineField1 const& function_1_coef_expected,
            SplineField2 const& function_2_coef_expected)
    {
        // --- get the spline representation on CPU
        auto function_coef_1_host = ddc::create_mirror_and_copy(std::get<0>(function_coef));
        auto function_coef_2_host = ddc::create_mirror_and_copy(std::get<1>(function_coef));

        // --- get expected spline representations on CPU
        auto function_coef_1_expected_host = ddc::create_mirror_and_copy(function_1_coef_expected);
        auto function_coef_2_expected_host = ddc::create_mirror_and_copy(function_2_coef_expected);

        // --- check error
        double max_abs_error = 0;
        ddc::for_each(get_idx_range(function_coef_1_expected_host), [&](auto const idx) {
            double err = abs(function_coef_1_expected_host(idx) - function_coef_1_host(idx));
            max_abs_error = std::max(max_abs_error, err);
        });
        ddc::for_each(get_idx_range(function_coef_2_expected_host), [&](auto const idx) {
            double err = abs(function_coef_2_expected_host(idx) - function_coef_2_host(idx));
            max_abs_error = std::max(max_abs_error, err);
        });
        EXPECT_LE(max_abs_error, 1e-15);
        std::cout << "Max absolute error: " << max_abs_error << std::endl;
    }
};


} // namespace


/*
    The global index range is splitted into two 1D patches on X1 and X2.
 */
TEST_F(MultipatchSplineBuilderTest, TwoPatches1D)
{
    // List of spline builders
    SplineX1Builder_1d builder_x1(idx_range_x1);
    SplineX2Builder_1d builder_x2(idx_range_x2);

    // Collection of builders for each patch
    MultipatchSplineBuilder builder(builder_x1, builder_x2);


    // Function tests
    // --- patch 1
    DFieldMem<IdxRangeX1> function_1_alloc(idx_range_x1);
    DField<IdxRangeX1> function_1 = get_field(function_1_alloc);

    // --- patch 2
    DFieldMem<IdxRangeX2> function_2_alloc(idx_range_x2);
    DField<IdxRangeX2> function_2 = get_field(function_2_alloc);

    initialize_1D_functions(function_1, function_2);

    // --- collection of values of the function on each patch
    std::tuple function_values = {function_1, function_2};



    // Spline representations
    // --- patch 1
    DFieldMem<IdxRange<BSplinesX1>> function_1_coef_alloc(get_spline_idx_range(builder_x1));
    DField<IdxRange<BSplinesX1>> function_1_coef = get_field(function_1_coef_alloc);

    DFieldMem<IdxRange<BSplinesX1>> function_1_coef_expected_alloc(
            get_spline_idx_range(builder_x1));
    DField<IdxRange<BSplinesX1>> function_1_coef_expected
            = get_field(function_1_coef_expected_alloc);

    // --- patch 2
    DFieldMem<IdxRange<BSplinesX2>> function_2_coef_alloc(get_spline_idx_range(builder_x2));
    DField<IdxRange<BSplinesX2>> function_2_coef = get_field(function_2_coef_alloc);

    DFieldMem<IdxRange<BSplinesX2>> function_2_coef_expected_alloc(
            get_spline_idx_range(builder_x2));
    DField<IdxRange<BSplinesX2>> function_2_coef_expected
            = get_field(function_2_coef_expected_alloc);

    // --- collection
    std::tuple function_coef = {function_1_coef, function_2_coef};


    // Build the spline representations on each patch
    builder(function_coef, function_values);



    // Check the results are the same than with one by one spline building
    builder_x1(function_1_coef_expected, get_const_field(function_1));
    builder_x2(function_2_coef_expected, get_const_field(function_2));
    check_if_equal_to_expected(function_coef, function_1_coef_expected, function_2_coef_expected);
}



/*
    The global index range is splitted into two 2D patches on (X1, Y1) and (X2, Y2).
 */
TEST_F(MultipatchSplineBuilderTest, TwoPatches2D)
{
    // List of index ranges
    IdxRangeXY1 idx_range_xy1(idx_range_x1, idx_range_y1);
    IdxRangeXY2 idx_range_xy2(idx_range_x2, idx_range_y2);


    // List of spline builders
    SplineX1Builder_2d builder_x1(idx_range_xy1);
    SplineX2Builder_2d builder_x2(idx_range_xy2);

    // Collection of builders for each patch
    MultipatchSplineBuilder builder(builder_x1, builder_x2);


    // Function tests
    // --- patch 1
    DFieldMem<IdxRangeXY1> function_1_alloc(idx_range_xy1);
    DField<IdxRangeXY1> function_1 = get_field(function_1_alloc);

    // --- patch 2
    DFieldMem<IdxRangeXY2> function_2_alloc(idx_range_xy2);
    DField<IdxRangeXY2> function_2 = get_field(function_2_alloc);

    initialize_2D_functions(function_1, function_2);

    // --- collection of values of the function on each patch
    std::tuple function_values = {function_1, function_2};



    // Spline representations
    // --- patch 1
    DFieldMem<IdxRange<BSplinesX1, GridY1>> function_1_coef_alloc(
            builder_x1.batched_spline_domain());
    DField<IdxRange<BSplinesX1, GridY1>> function_1_coef = get_field(function_1_coef_alloc);

    DFieldMem<IdxRange<BSplinesX1, GridY1>> function_1_coef_expected_alloc(
            builder_x1.batched_spline_domain());
    DField<IdxRange<BSplinesX1, GridY1>> function_1_coef_expected
            = get_field(function_1_coef_expected_alloc);

    // --- patch 2
    DFieldMem<IdxRange<BSplinesX2, GridY2>> function_2_coef_alloc(
            builder_x2.batched_spline_domain());
    DField<IdxRange<BSplinesX2, GridY2>> function_2_coef = get_field(function_2_coef_alloc);

    DFieldMem<IdxRange<BSplinesX2, GridY2>> function_2_coef_expected_alloc(
            builder_x2.batched_spline_domain());
    DField<IdxRange<BSplinesX2, GridY2>> function_2_coef_expected
            = get_field(function_2_coef_expected_alloc);

    // --- collection
    std::tuple function_coef = {function_1_coef, function_2_coef};


    // Build the spline representations on each patch
    builder(function_coef, function_values);



    // Check the results are the same than with one by one spline building
    builder_x1(function_1_coef_expected, get_const_field(function_1));
    builder_x2(function_2_coef_expected, get_const_field(function_2));
    check_if_equal_to_expected(function_coef, function_1_coef_expected, function_2_coef_expected);
}
