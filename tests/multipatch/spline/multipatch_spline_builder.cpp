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
struct RDimX1
{
    static bool constexpr PERIODIC = false;
};

struct RDimY1
{
    static bool constexpr PERIODIC = false;
};

// Continuous dimension of patch 2
struct RDimX2
{
    static bool constexpr PERIODIC = false;
};

struct RDimY2
{
    static bool constexpr PERIODIC = false;
};



// Splines
struct BSplinesX1 : ddc::UniformBSplines<RDimX1, 3>
{
};
struct BSplinesY1 : ddc::UniformBSplines<RDimY1, 3>
{
};
struct BSplinesX2 : ddc::UniformBSplines<RDimX2, 3>
{
};
struct BSplinesY2 : ddc::UniformBSplines<RDimY2, 3>
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
struct IDimX1 : SplineInterpPointsX1::interpolation_mesh_type
{
};
struct IDimY1 : SplineInterpPointsY1::interpolation_mesh_type
{
};
// Discrete dimension of patch 2
struct IDimX2 : SplineInterpPointsX2::interpolation_mesh_type
{
};
struct IDimY2 : SplineInterpPointsY2::interpolation_mesh_type
{
};



/*
    TODO: use Patch class once merged. 
*/
using CoordX1 = ddc::Coordinate<RDimX1>;
using CoordY1 = ddc::Coordinate<RDimY1>;
using CoordX2 = ddc::Coordinate<RDimX2>;
using CoordY2 = ddc::Coordinate<RDimY2>;


using IndexX1 = ddc::DiscreteElement<IDimX1>;
using IndexY1 = ddc::DiscreteElement<IDimY1>;
using IndexX2 = ddc::DiscreteElement<IDimX2>;
using IndexY2 = ddc::DiscreteElement<IDimY2>;

using IndexXY1 = ddc::DiscreteElement<IDimX1, IDimY1>;
using IndexXY2 = ddc::DiscreteElement<IDimX2, IDimY2>;


using IVectX1 = ddc::DiscreteVector<IDimX1>;
using IVectY1 = ddc::DiscreteVector<IDimY1>;
using IVectX2 = ddc::DiscreteVector<IDimX2>;
using IVectY2 = ddc::DiscreteVector<IDimY2>;


using IDomainX1 = ddc::DiscreteDomain<IDimX1>;
using IDomainY1 = ddc::DiscreteDomain<IDimY1>;
using IDomainX2 = ddc::DiscreteDomain<IDimX2>;
using IDomainY2 = ddc::DiscreteDomain<IDimY2>;

using IDomainXY1 = ddc::DiscreteDomain<IDimX1, IDimY1>;
using IDomainXY2 = ddc::DiscreteDomain<IDimX2, IDimY2>;


// Operators
using SplineX1Builder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX1,
        IDimX1,
        SplineX1Boundary,
        SplineX1Boundary,
        ddc::SplineSolver::GINKGO,
        IDimX1>;

using SplineX2Builder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX2,
        IDimX2,
        SplineX2Boundary,
        SplineX2Boundary,
        ddc::SplineSolver::GINKGO,
        IDimX2>;



using SplineX1Builder_2d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX1,
        IDimX1,
        SplineX1Boundary,
        SplineX1Boundary,
        ddc::SplineSolver::GINKGO,
        IDimX1,
        IDimY1>;

using SplineX2Builder_2d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX2,
        IDimX2,
        SplineX2Boundary,
        SplineX2Boundary,
        ddc::SplineSolver::GINKGO,
        IDimX2,
        IDimY2>;



class MultipatchSplineBuilderTest : public ::testing::Test
{
private:
    static constexpr IVectX1 x1_size = IVectX1(10);
    static constexpr IVectY1 y1_size = IVectY1(10);

    static constexpr IVectX2 x2_size = IVectX2(10);
    static constexpr IVectY2 y2_size = IVectY2(10);

protected:
    IDomainX1 const domain_x1;
    IDomainY1 const domain_y1;

    IDomainX2 const domain_x2;
    IDomainY2 const domain_y2;

public:
    MultipatchSplineBuilderTest()
        : domain_x1(SplineInterpPointsX1::get_domain<IDimX1>())
        , domain_y1(SplineInterpPointsY1::get_domain<IDimY1>())
        , domain_x2(SplineInterpPointsX2::get_domain<IDimX2>())
        , domain_y2(SplineInterpPointsY2::get_domain<IDimY2>()) {};

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

        ddc::init_discrete_space<IDimX1>(SplineInterpPointsX1::get_sampling<IDimX1>());
        ddc::init_discrete_space<IDimY1>(SplineInterpPointsY1::get_sampling<IDimY1>());

        // Patch 2
        CoordX2 const x2_min(2 * M_PI);
        CoordX2 const x2_max(4 * M_PI);

        CoordY2 const y2_min(0.0);
        CoordY2 const y2_max(1.0);

        ddc::init_discrete_space<BSplinesX2>(x2_min, x2_max, x2_size);
        ddc::init_discrete_space<BSplinesY2>(y2_min, y2_max, y2_size);

        ddc::init_discrete_space<IDimX2>(SplineInterpPointsX2::get_sampling<IDimX2>());
        ddc::init_discrete_space<IDimY2>(SplineInterpPointsY2::get_sampling<IDimY2>());
    }


    // Initialisation methods --------------------------------------------------------------------
    void initialize_1D_functions(
            device_t<ddc::ChunkSpan<double, IDomainX1>> function_1,
            device_t<ddc::ChunkSpan<double, IDomainX2>> function_2)
    {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                domain_x1,
                KOKKOS_LAMBDA(IndexX1 const idx) {
                    double const x = ddc::coordinate(idx);
                    function_1(idx) = Kokkos::sin(x);
                });

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                domain_x2,
                KOKKOS_LAMBDA(IndexX2 const idx) {
                    double const x = ddc::coordinate(idx);
                    function_2(idx) = Kokkos::sin(2 * x);
                });
    }


    void initialize_2D_functions(
            device_t<ddc::ChunkSpan<double, IDomainXY1>> function_1,
            device_t<ddc::ChunkSpan<double, IDomainXY2>> function_2)
    {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                function_1.domain(),
                KOKKOS_LAMBDA(IndexXY1 const idx) {
                    double const x = ddc::coordinate(IndexX1(idx));
                    double const y = ddc::coordinate(IndexY1(idx));
                    function_1(idx) = Kokkos::sin(x) + Kokkos::cos(y);
                });

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                function_2.domain(),
                KOKKOS_LAMBDA(IndexXY2 const idx) {
                    double const x = ddc::coordinate(IndexX2(idx));
                    double const y = ddc::coordinate(IndexY2(idx));
                    function_2(idx) = Kokkos::sin(2 * x) + Kokkos::cos(y);
                });
    }


    // Test method -------------------------------------------------------------------------------
    template <class CoefTuple, class SplineSpan1, class SplineSpan2>
    void check_if_equal_to_expected(
            CoefTuple const& function_coef,
            SplineSpan1 const& function_1_coef_expected,
            SplineSpan2 const& function_2_coef_expected)
    {
        // --- get the spline representation on CPU
        auto function_coef_1_host = ddc::create_mirror_and_copy(std::get<0>(function_coef));
        auto function_coef_2_host = ddc::create_mirror_and_copy(std::get<1>(function_coef));

        // --- get expected spline representations on CPU
        auto function_coef_1_expected_host = ddc::create_mirror_and_copy(function_1_coef_expected);
        auto function_coef_2_expected_host = ddc::create_mirror_and_copy(function_2_coef_expected);

        // --- check error
        double max_abs_error = 0;
        ddc::for_each(function_coef_1_expected_host.domain(), [&](auto const idx) {
            double err = abs(function_coef_1_expected_host(idx) - function_coef_1_host(idx));
            max_abs_error = std::max(max_abs_error, err);
        });
        ddc::for_each(function_coef_2_expected_host.domain(), [&](auto const idx) {
            double err = abs(function_coef_2_expected_host(idx) - function_coef_2_host(idx));
            max_abs_error = std::max(max_abs_error, err);
        });
        EXPECT_LE(max_abs_error, 1e-15);
        std::cout << "Max absolute error: " << max_abs_error << std::endl;
    }
};


} // namespace


/*
    The global domain is splitted into two 1D patches on RDimX1 and RDimX2.
 */
TEST_F(MultipatchSplineBuilderTest, TwoPatches1D)
{
    // List of spline builders
    SplineX1Builder_1d builder_x1(domain_x1);
    SplineX2Builder_1d builder_x2(domain_x2);

    // Collection of builders for each patch
    MultipatchSplineBuilder builder(builder_x1, builder_x2);


    // Function tests
    // --- patch 1
    device_t<ddc::Chunk<double, IDomainX1>> function_1_alloc(domain_x1);
    ddc::ChunkSpan function_1 = function_1_alloc.span_view();

    // --- patch 2
    device_t<ddc::Chunk<double, IDomainX2>> function_2_alloc(domain_x2);
    ddc::ChunkSpan function_2 = function_2_alloc.span_view();

    initialize_1D_functions(function_1, function_2);

    // --- collection of values of the function on each patch
    std::tuple function_values = {function_1, function_2};



    // Spline representations
    // --- patch 1
    device_t<ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX1>>> function_1_coef_alloc(
            builder_x1.spline_domain());
    ddc::ChunkSpan function_1_coef = function_1_coef_alloc.span_view();

    device_t<ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX1>>> function_1_coef_expected_alloc(
            builder_x1.spline_domain());
    ddc::ChunkSpan function_1_coef_expected = function_1_coef_expected_alloc.span_view();

    // --- patch 2
    device_t<ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX2>>> function_2_coef_alloc(
            builder_x2.spline_domain());
    ddc::ChunkSpan function_2_coef = function_2_coef_alloc.span_view();

    device_t<ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX2>>> function_2_coef_expected_alloc(
            builder_x2.spline_domain());
    ddc::ChunkSpan function_2_coef_expected = function_2_coef_expected_alloc.span_view();

    // --- collection
    std::tuple function_coef = {function_1_coef, function_2_coef};


    // Build the spline representations on each patch
    builder(function_coef, function_values);



    // Check the results are the same than with one by one spline building
    builder_x1(function_1_coef_expected, function_1.span_cview());
    builder_x2(function_2_coef_expected, function_2.span_cview());
    check_if_equal_to_expected(function_coef, function_1_coef_expected, function_2_coef_expected);
}



/*
    The global domain is splitted into two 2D patches on (RDimX1, RDimY1) and (RDimX2, RDimY2).
 */
TEST_F(MultipatchSplineBuilderTest, TwoPatches2D)
{
    // List of domains
    IDomainXY1 domain_xy1(domain_x1, domain_y1);
    IDomainXY2 domain_xy2(domain_x2, domain_y2);


    // List of spline builders
    SplineX1Builder_2d builder_x1(domain_xy1);
    SplineX2Builder_2d builder_x2(domain_xy2);

    // Collection of builders for each patch
    MultipatchSplineBuilder builder(builder_x1, builder_x2);


    // Function tests
    // --- patch 1
    device_t<ddc::Chunk<double, IDomainXY1>> function_1_alloc(domain_xy1);
    ddc::ChunkSpan function_1 = function_1_alloc.span_view();

    // --- patch 2
    device_t<ddc::Chunk<double, IDomainXY2>> function_2_alloc(domain_xy2);
    ddc::ChunkSpan function_2 = function_2_alloc.span_view();

    initialize_2D_functions(function_1, function_2);

    // --- collection of values of the function on each patch
    std::tuple function_values = {function_1, function_2};



    // Spline representations
    // --- patch 1
    device_t<ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX1, IDimY1>>> function_1_coef_alloc(
            builder_x1.spline_domain());
    ddc::ChunkSpan function_1_coef = function_1_coef_alloc.span_view();

    device_t<ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX1, IDimY1>>>
            function_1_coef_expected_alloc(builder_x1.spline_domain());
    ddc::ChunkSpan function_1_coef_expected = function_1_coef_expected_alloc.span_view();

    // --- patch 2
    device_t<ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX2, IDimY2>>> function_2_coef_alloc(
            builder_x2.spline_domain());
    ddc::ChunkSpan function_2_coef = function_2_coef_alloc.span_view();

    device_t<ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX2, IDimY2>>>
            function_2_coef_expected_alloc(builder_x2.spline_domain());
    ddc::ChunkSpan function_2_coef_expected = function_2_coef_expected_alloc.span_view();

    // --- collection
    std::tuple function_coef = {function_1_coef, function_2_coef};


    // Build the spline representations on each patch
    builder(function_coef, function_values);



    // Check the results are the same than with one by one spline building
    builder_x1(function_1_coef_expected, function_1.span_cview());
    builder_x2(function_2_coef_expected, function_2.span_cview());
    check_if_equal_to_expected(function_coef, function_1_coef_expected, function_2_coef_expected);
}
