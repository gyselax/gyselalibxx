// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_non_uniform.hpp"
#include "connectivity.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "mesh_builder.hpp"
#include "multipatch_field.hpp"
#include "multipatch_field_mem.hpp"
#include "spline_builder_2d_collection.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"


namespace {
using namespace onion_shape_non_uniform_2d_2patches;

// INTERFACES ------------------------------------------------------------------------------------
using NorthInterface1 = Interface<NorthEdge<1>, OutsideEdge, true>;
using SouthInterface1 = Interface<SouthEdge<1>, OutsideEdge, true>;
using WestInterface1 = Interface<WestEdge<1>, OutsideEdge, true>;

using NorthInterface2 = Interface<NorthEdge<2>, OutsideEdge, true>;
using SouthInterface2 = Interface<SouthEdge<2>, OutsideEdge, true>;
using EastInterface2 = Interface<EastEdge<2>, OutsideEdge, true>;

using Interface_1_2 = Interface<EastEdge<1>, WestEdge<2>, true>;

// CONNECTIVITY ----------------------------------------------------------------------------------
using Connectivity = MultipatchConnectivity<
        NorthInterface1,
        SouthInterface1,
        WestInterface1,
        NorthInterface2,
        SouthInterface2,
        EastInterface2,
        Interface_1_2>;

template <int PatchIdx>
using SplineInterpPointsR = ddc::GrevilleInterpolationPoints<
        BSplinesR<PatchIdx>,
        ddc::SplineBuilderClosure::HERMITE,
        ddc::SplineBuilderClosure::HERMITE>;
template <int PatchIdx>
using SplineInterpPointsTheta = ddc::KnotsAsInterpolationPoints<
        BSplinesTheta<PatchIdx>,
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineBuilderClosure::PERIODIC>;

// Operators
template <int PatchIdx>
using SplineBuilderOnPatch = ddc::SplineBuilder2D<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesR<PatchIdx>,
        BSplinesTheta<PatchIdx>,
        GridR<PatchIdx>,
        GridTheta<PatchIdx>,
        ddc::SplineBuilderClosure::HERMITE,
        ddc::SplineBuilderClosure::HERMITE,
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineSolver::LAPACK>;

using SplineBuilderCollection = SplineBuilder2DCollection<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplines1OnPatch,
        BSplines2OnPatch,
        Grid1OnPatch,
        Grid2OnPatch,
        ddc::SplineBuilderClosure::HERMITE,
        ddc::SplineBuilderClosure::HERMITE,
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineBuilderClosure::HERMITE,
        Connectivity,
        ddc::SplineSolver::LAPACK,
        DConstFieldOnPatch,
        Patch1,
        Patch2>;

class SplineBuilderCollection2DOnionTest : public ::testing::Test
{
public:
    SplineBuilderCollection2DOnionTest() {}

    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        Coord<R> const r1_min(0.0);
        Coord<R> const r1_max(1.0);
        IdxStep<GridR<1>> r1_ncells(10);

        Coord<Theta> const theta1_min(0.0);
        Coord<Theta> const theta1_max(2 * M_PI);
        IdxStep<GridTheta<1>> theta1_ncells(12);

        ddc::init_discrete_space<BSplinesR<1>>(
                build_uniform_break_points(r1_min, r1_max, r1_ncells));
        ddc::init_discrete_space<BSplinesTheta<1>>(
                build_uniform_break_points(theta1_min, theta1_max, theta1_ncells));

        // Patch 2
        Coord<R> const r2_min(1.0);
        Coord<R> const r2_max(2.0);
        IdxStep<GridR<2>> r2_ncells(6);

        Coord<Theta> const theta2_min(0.0);
        Coord<Theta> const theta2_max(2 * M_PI);
        IdxStep<GridTheta<2>> theta2_ncells(8);

        ddc::init_discrete_space<BSplinesR<2>>(
                build_uniform_break_points(r2_min, r2_max, r2_ncells));
        ddc::init_discrete_space<BSplinesTheta<2>>(
                build_uniform_break_points(theta2_min, theta2_max, theta2_ncells));
    }

    template <class ExecSpace, class Layout>
    void initialise_2D_functions(
            DField<typename Patch1::IdxRange12, ExecSpace, Layout> function_1,
            DField<typename Patch2::IdxRange12, ExecSpace, Layout> function_2)
    {
        ddc::parallel_for_each(
                get_idx_range(function_1),
                KOKKOS_LAMBDA(Idx<GridR<1>, GridTheta<1>> idx) {
                    double const r = ddc::coordinate(Idx<GridR<1>>(idx));
                    double const theta = ddc::coordinate(Idx<GridTheta<1>>(idx));
                    function_1(idx) = Kokkos::sin(r) + Kokkos::cos(theta);
                });
        ddc::parallel_for_each(
                get_idx_range(function_2),
                KOKKOS_LAMBDA(Idx<GridR<2>, GridTheta<2>> idx) {
                    double const r = ddc::coordinate(Idx<GridR<2>>(idx));
                    double const theta = ddc::coordinate(Idx<GridTheta<2>>(idx));
                    function_2(idx) = Kokkos::sin(2 * r * theta);
                });
    }

    template <class ExecSpace, class Layout>
    void initialise_2D_derivatives_2(
            DField<IdxRange<ddc::Deriv<R>, GridTheta<2>>, ExecSpace, Layout> derivs_rmin2,
            DField<IdxRange<ddc::Deriv<R>, GridTheta<2>>, ExecSpace, Layout> derivs_rmax2)
    {
        double const rmin2 = ddc::discrete_space<typename Patch2::BSplines1>().rmin();
        double const rmax2 = ddc::discrete_space<typename Patch2::BSplines1>().rmax();
        ddc::parallel_for_each(
                get_idx_range(derivs_rmin2),
                KOKKOS_LAMBDA(Idx<ddc::Deriv<R>, GridTheta<2>> idx) {
                    double const theta = ddc::coordinate(Idx<GridTheta<2>>(idx));
                    derivs_rmin2(idx) = 2 * theta * Kokkos::cos(2 * rmin2 * theta);
                    derivs_rmax2(idx) = 2 * theta * Kokkos::cos(2 * rmax2 * theta);
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
        ddc::host_for_each(
                get_idx_range(function_coef_1_expected_host),
                [&](Idx<Patch1::BSplines1, Patch1::BSplines2> const idx) {
                    double err = std::abs(
                            function_coef_1_expected_host(idx) - function_coef_1_host(idx));
                    max_abs_error = std::max(max_abs_error, err);
                });
        ddc::host_for_each(
                get_idx_range(function_coef_2_expected_host),
                [&](Idx<Patch2::BSplines1, Patch2::BSplines2> const idx) {
                    double err = std::abs(
                            function_coef_2_expected_host(idx) - function_coef_2_host(idx));
                    max_abs_error = std::max(max_abs_error, err);
                });
        std::cout << "Max absolute error: " << max_abs_error << std::endl;
        EXPECT_LE(max_abs_error, 1e-12);
    }
};

} // namespace


TEST_F(SplineBuilderCollection2DOnionTest, SameContinuousDimTest)
{
    ddc::init_discrete_space<GridR<1>>(SplineInterpPointsR<1>::get_sampling<GridR<1>>());
    ddc::init_discrete_space<GridTheta<1>>(
            SplineInterpPointsTheta<1>::get_sampling<GridTheta<1>>());

    ddc::init_discrete_space<GridR<2>>(SplineInterpPointsR<2>::get_sampling<GridR<2>>());
    ddc::init_discrete_space<GridTheta<2>>(
            SplineInterpPointsTheta<2>::get_sampling<GridTheta<2>>());

    IdxRange<GridR<1>> const idx_range_r1 = SplineInterpPointsR<1>::get_domain<GridR<1>>();
    IdxRange<GridTheta<1>> const idx_range_theta1
            = SplineInterpPointsTheta<1>::get_domain<GridTheta<1>>();

    IdxRange<GridR<2>> const idx_range_r2 = SplineInterpPointsR<2>::get_domain<GridR<2>>();
    IdxRange<GridTheta<2>> const idx_range_theta2
            = SplineInterpPointsTheta<2>::get_domain<GridTheta<2>>();

    IdxRange<GridR<1>, GridTheta<1>> const idx_range_rtheta1(idx_range_r1, idx_range_theta1);
    IdxRange<GridR<2>, GridTheta<2>> const idx_range_rtheta2(idx_range_r2, idx_range_theta2);

    // List of spline builders
    SplineBuilderOnPatch<1> builder1(idx_range_rtheta1);
    SplineBuilderOnPatch<2> builder2(idx_range_rtheta2);

    // Collection of builders for each patch
    SplineBuilderCollection builders(builder1, builder2);

    // Function tests
    // --- patch 1
    DFieldMemOnPatch<Patch1> function_1_alloc(idx_range_rtheta1);
    DFieldOnPatch<Patch1> function_1 = get_field(function_1_alloc);

    // --- patch 2
    DFieldMemOnPatch<Patch2> function_2_alloc(idx_range_rtheta2);
    DFieldOnPatch<Patch2> function_2 = get_field(function_2_alloc);

    initialise_2D_functions(function_1, function_2);

    // --- collection of values of the function on each patch
    MultipatchField<DConstFieldOnPatch, Patch1, Patch2>
            function_values(get_const_field(function_1), get_const_field(function_2));

    // Derivatives Patch1
    IdxRange<ddc::Deriv<R>> deriv_r_idx_range(Idx<ddc::Deriv<R>>(1), IdxStep<ddc::Deriv<R>>(1));

    IdxRange<ddc::Deriv<R>, GridTheta<1>> derivs_r1_idx_range(deriv_r_idx_range, idx_range_theta1);
    DFieldMem<IdxRange<ddc::Deriv<R>, GridTheta<1>>> derivs_rmin1_alloc(derivs_r1_idx_range);
    DFieldMem<IdxRange<ddc::Deriv<R>, GridTheta<1>>> derivs_rmax1_alloc(derivs_r1_idx_range);
    DField<IdxRange<ddc::Deriv<R>, GridTheta<1>>> derivs_rmin1 = get_field(derivs_rmin1_alloc);
    DField<IdxRange<ddc::Deriv<R>, GridTheta<1>>> derivs_rmax1 = get_field(derivs_rmax1_alloc);

    double const rmin1 = ddc::coordinate(idx_range_r1.front());
    double const rmax1 = ddc::coordinate(idx_range_r1.back());
    ddc::parallel_fill(derivs_rmin1, Kokkos::cos(rmin1));
    ddc::parallel_fill(derivs_rmax1, Kokkos::cos(rmax1));

    // Derivatives Patch 2
    IdxRange<ddc::Deriv<R>, GridTheta<2>> derivs_r2_idx_range(deriv_r_idx_range, idx_range_theta2);
    DFieldMem<IdxRange<ddc::Deriv<R>, GridTheta<2>>> derivs_rmin2_alloc(derivs_r2_idx_range);
    DFieldMem<IdxRange<ddc::Deriv<R>, GridTheta<2>>> derivs_rmax2_alloc(derivs_r2_idx_range);
    DField<IdxRange<ddc::Deriv<R>, GridTheta<2>>> derivs_rmin2 = get_field(derivs_rmin2_alloc);
    DField<IdxRange<ddc::Deriv<R>, GridTheta<2>>> derivs_rmax2 = get_field(derivs_rmax2_alloc);

    initialise_2D_derivatives_2(derivs_rmin2, derivs_rmax2);

    // --- collection of values of the function on each patch
    MultipatchField<ConstDeriv1_OnPatch_2D, Patch1, Patch2>
            derivs_rmin(get_const_field(derivs_rmin1), get_const_field(derivs_rmin2));
    MultipatchField<ConstDeriv1_OnPatch_2D, Patch1, Patch2>
            derivs_rmax(get_const_field(derivs_rmax1), get_const_field(derivs_rmax2));

    // Spline representations
    // --- patch 1
    DFieldMem<typename Patch1::IdxRangeBS12> function_1_coef_expected_alloc(
            builder1.batched_spline_domain(idx_range_rtheta1));
    SplineCoeffOnPatch_2D<Patch1> function_1_coef_expected
            = get_field(function_1_coef_expected_alloc);

    // --- patch 2
    DFieldMem<typename Patch2::IdxRangeBS12> function_2_coef_expected_alloc(
            builder2.batched_spline_domain(idx_range_rtheta2));
    SplineCoeffOnPatch_2D<Patch2> function_2_coef_expected
            = get_field(function_2_coef_expected_alloc);

    // --- collection
    MultipatchFieldMem<SplineCoeffMemOnPatch_2D, Patch1, Patch2> function_coef_alloc(
            builders.spline_idx_ranges());
    MultipatchField<SplineCoeffOnPatch_2D, Patch1, Patch2> function_coef(function_coef_alloc);

    // Build the spline representations on each patch
    builders(function_coef, function_values, derivs_rmin, derivs_rmax);

    // Check the results are the same than with one by one spline building
    builder1(
            function_1_coef_expected,
            get_const_field(function_1),
            std::optional(get_const_field(derivs_rmin1)),
            std::optional(get_const_field(derivs_rmax1)));
    builder2(
            function_2_coef_expected,
            get_const_field(function_2),
            std::optional(get_const_field(derivs_rmin2)),
            std::optional(get_const_field(derivs_rmax2)));

    check_if_equal_to_expected(function_coef, function_1_coef_expected, function_2_coef_expected);
}
