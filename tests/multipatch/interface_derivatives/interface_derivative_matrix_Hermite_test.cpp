// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "3patches_2d_non_periodic_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "interface.hpp"
#include "interface_derivative_matrix.hpp"
#include "interface_derivatives_matrix_test_utils.hpp"
#include "interface_derivatives_test_utils.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "single_interface_derivatives_calculator.hpp"
#include "single_interface_derivatives_calculator_collection.hpp"
#include "types.hpp"
#include "view.hpp"


/*
    Test InterfaceDerivativeMatrix on the following geometry:

        |  1  |  2  |  3  | 

        with the global X dimension with Hermite boundary conditions 
        and the global Y spline with additional points as closure condition 
        (ddc::BoundCond::GREVILLE).

    > test ddc::BoundCond::HERMITE boundary conditions. 
    > test application on the X direction. 
    > test application to compute first derivatives and cross-derivatives. 
    > test on a non-uniform patches. 
    > test with exact formulation in SingleInterfaceDerivativeCalculator. 
    > test agreement between computed and global spline derivatives. 
    > test agreement between local and global splines.
*/

namespace {
// Multi-patch tags ---
// using namespace periodic_strips_non_uniform_2d_9patches;
using namespace non_periodic_non_uniform_2d_3patches;

// INTERFACES ------------------------------------------------------------------------------------
using NorthInterface1 = Interface<NorthEdge<1>, OutsideEdge, true>;
using NorthInterface2 = Interface<NorthEdge<2>, OutsideEdge, true>;
using NorthInterface3 = Interface<NorthEdge<3>, OutsideEdge, true>;

using SouthInterface1 = Interface<OutsideEdge, SouthEdge<1>, true>;
using SouthInterface2 = Interface<OutsideEdge, SouthEdge<2>, true>;
using SouthInterface3 = Interface<OutsideEdge, SouthEdge<3>, true>;

using WestInterface1 = Interface<OutsideEdge, WestEdge<1>, true>;
using EastInterface3 = Interface<EastEdge<3>, OutsideEdge, true>;

using Interface_1_2 = Interface<EastEdge<1>, WestEdge<2>, true>;
using Interface_2_3 = Interface<EastEdge<2>, WestEdge<3>, true>;


// CONNECTIVITY ----------------------------------------------------------------------------------
using Connectivity = MultipatchConnectivity<
        NorthInterface1,
        NorthInterface2,
        NorthInterface3,
        SouthInterface1,
        SouthInterface2,
        SouthInterface3,
        WestInterface1,
        EastInterface3,
        Interface_1_2,
        Interface_2_3>;

// TODO: fix error here.
// using interface_collection = typename Connectivity::get_all_interfaces_along_direction_t<GridX<1>>;
// using all_patches = typename Connectivity::all_patches;

// static_assert(std::is_same_v<
//               interface_collection,
//               ddc::detail::TypeSeq<WestInterface1, Interface_1_2, Interface_2_3, EastInterface3>>);

// using Grid1DSeq = collect_grids_on_dim_t<Patch1, GridX<1>, interface_collection>;

// static_assert(std::is_same_v<Grid1DSeq, ddc::detail::TypeSeq<GridX<1>, GridX<2>, GridX<3>>>);
// static_assert(ddc::in_tags_v<GridX<1>, Grid1DSeq>);
// static_assert(ddc::in_tags_v<GridX<2>, Grid1DSeq>);
// static_assert(ddc::in_tags_v<GridX<3>, Grid1DSeq>);


// using Grid1DSeq2 = collect_grids_on_dim_t<Patch3, GridX<3>, interface_collection>;

// static_assert(std::is_same_v<Grid1DSeq2, ddc::detail::TypeSeq<GridX<1>, GridX<2>, GridX<3>>>);
// static_assert(ddc::in_tags_v<GridX<1>, Grid1DSeq2>);
// static_assert(ddc::in_tags_v<GridX<2>, Grid1DSeq2>);
// static_assert(ddc::in_tags_v<GridX<3>, Grid1DSeq2>);
// ---
// template <std::size_t Index>
// using PatchI = ddc::type_seq_element_t<Index, PatchOrdering>;

// Equivalent global mesh tags ---
struct Xg
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Xg;
};

struct Yg
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Yg;
};

struct GridXg : NonUniformGridBase<Xg>
{
};
struct GridYg : NonUniformGridBase<Yg>
{
};

struct BSplinesXg : ddc::NonUniformBSplines<Xg, 3>
{
};
struct BSplinesYg : ddc::NonUniformBSplines<Yg, 3>
{
};

using DerivXg = ddc::Deriv<Xg>;
using DerivYg = ddc::Deriv<Yg>;

using HostExecSpace = Kokkos::DefaultHostExecutionSpace;


// Interpolation points type for the patches.
template <std::size_t PatchIdx>
using SplineInterpPointsX = ddcHelper::NonUniformInterpolationPoints<
        BSplinesX<PatchIdx>,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

template <std::size_t PatchIdx>
using SplineInterpPointsY = ddcHelper::NonUniformInterpolationPoints<
        BSplinesY<PatchIdx>,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

// Interpolation points type for the equivalent global spline.
using SplineInterpPointsXg = ddcHelper::
        NonUniformInterpolationPoints<BSplinesXg, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;
using SplineInterpPointsYg = ddcHelper::
        NonUniformInterpolationPoints<BSplinesYg, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;


// Operators on the equivalent global spline.
using SplineRThetagBuilder = ddc::SplineBuilder2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesXg,
        BSplinesYg,
        GridXg,
        GridYg,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::SplineSolver::LAPACK>;

using SplineRThetagEvaluator = ddc::SplineEvaluator2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesXg,
        BSplinesYg,
        GridXg,
        GridYg,
        ddc::ConstantExtrapolationRule<Xg, Yg>,
        ddc::ConstantExtrapolationRule<Xg, Yg>,
        ddc::ConstantExtrapolationRule<Yg, Xg>,
        ddc::ConstantExtrapolationRule<Yg, Xg>>;



struct InterfaceDerivativeMatrixHermiteTest : public ::testing::Test
{
    // DEFINE BOUNDARIES OF THE DOMAINS ----------------------------------------------------------
    // // patches 1 ---------------------------------
    static constexpr Coord<X<1>> x1_min = Coord<X<1>>(0.0);
    static constexpr Coord<X<1>> x1_max = Coord<X<1>>(1.0);
    static constexpr IdxStep<GridX<1>> x1_ncells = IdxStep<GridX<1>>(6);

    static constexpr Coord<Y<1>> y1_min = Coord<Y<1>>(0.0);
    static constexpr Coord<Y<1>> y1_max = Coord<Y<1>>(1.0);
    static constexpr IdxStep<GridY<1>> y1_ncells = IdxStep<GridY<1>>(6);

    // patches 2 ---------------------------------
    static constexpr Coord<X<2>> x2_min = Coord<X<2>>(1.0);
    static constexpr Coord<X<2>> x2_max = Coord<X<2>>(2.0);
    static constexpr IdxStep<GridX<2>> x2_ncells = IdxStep<GridX<2>>(6);

    static constexpr Coord<Y<2>> y2_min = Coord<Y<2>>(0.0);
    static constexpr Coord<Y<2>> y2_max = Coord<Y<2>>(1.0);
    static constexpr IdxStep<GridY<2>> y2_ncells = IdxStep<GridY<2>>(6);

    // patches 3 ---------------------------------
    static constexpr Coord<X<3>> x3_min = Coord<X<3>>(2.0);
    static constexpr Coord<X<3>> x3_max = Coord<X<3>>(3.0);
    static constexpr IdxStep<GridX<3>> x3_ncells = IdxStep<GridX<3>>(6);

    static constexpr Coord<Y<3>> y3_min = Coord<Y<3>>(0.0);
    static constexpr Coord<Y<3>> y3_max = Coord<Y<3>>(1.0);
    static constexpr IdxStep<GridY<3>> y3_ncells = IdxStep<GridY<3>>(6);

    // global ------------------------------------
    static constexpr Coord<Xg> xg_min = Coord<Xg> {double(x1_min)};
    static constexpr Coord<Xg> xg_max = Coord<Xg> {double(x3_max)};
    static constexpr IdxStep<GridXg> xg_ncells
            = IdxStep<GridXg>(x1_ncells.value() + x2_ncells.value() + x3_ncells.value());

    static constexpr Coord<Yg> yg_min = Coord<Yg> {double(y1_min)};
    static constexpr Coord<Yg> yg_max = Coord<Yg> {double(y1_max)};
    static constexpr IdxStep<GridYg> yg_ncells = IdxStep<GridYg>(y1_ncells.value());

protected:
    const IdxRange<GridX<1>> idx_range_x1;
    const IdxRange<GridX<2>> idx_range_x2;
    const IdxRange<GridX<3>> idx_range_x3;

    const IdxRange<GridY<1>> idx_range_y1;
    const IdxRange<GridY<2>> idx_range_y2;
    const IdxRange<GridY<3>> idx_range_y3;

    const IdxRange<GridX<1>, GridY<1>> idx_range_xy1;
    const IdxRange<GridX<2>, GridY<2>> idx_range_xy2;
    const IdxRange<GridX<3>, GridY<3>> idx_range_xy3;

    const IdxRange<GridXg> idx_range_xg;
    const IdxRange<GridYg> idx_range_yg;
    const IdxRange<GridXg, GridYg> idx_range_xy_g;

public:
    InterfaceDerivativeMatrixHermiteTest()
        : idx_range_x1(SplineInterpPointsX<1>::template get_domain<GridX<1>>())
        , idx_range_x2(SplineInterpPointsX<2>::template get_domain<GridX<2>>())
        , idx_range_x3(SplineInterpPointsX<3>::template get_domain<GridX<3>>())
        , idx_range_y1(SplineInterpPointsY<1>::template get_domain<GridY<1>>())
        , idx_range_y2(SplineInterpPointsY<2>::template get_domain<GridY<2>>())
        , idx_range_y3(SplineInterpPointsY<3>::template get_domain<GridY<3>>())
        , idx_range_xy1(idx_range_x1, idx_range_y1)
        , idx_range_xy2(idx_range_x2, idx_range_y2)
        , idx_range_xy3(idx_range_x3, idx_range_y3)
        , idx_range_xg(SplineInterpPointsXg::get_domain<GridXg>())
        , idx_range_yg(SplineInterpPointsYg::get_domain<GridYg>())
        , idx_range_xy_g(idx_range_xg, idx_range_yg)
    {
    }

    // INITIALISE DOMAINS ------------------------------------------------------------------------
    static void SetUpTestSuite()
    {
        // Creating of meshes and supports .......................................................
        // The patches are conforming between each other.
        std::vector<Coord<X<1>>> break_points_x1
                = build_random_non_uniform_break_points(x1_min, x1_max, x1_ncells);
        std::vector<Coord<X<2>>> break_points_x2
                = build_random_non_uniform_break_points(x2_min, x2_max, x2_ncells);
        std::vector<Coord<X<3>>> break_points_x3
                = build_random_non_uniform_break_points(x3_min, x3_max, x3_ncells);

        std::vector<Coord<Y<1>>> break_points_y1
                = build_random_non_uniform_break_points(y1_min, y1_max, y1_ncells);
        std::vector<Coord<Y<2>>> break_points_y2
                = build_random_non_uniform_break_points(y2_min, y2_max, y2_ncells);
        std::vector<Coord<Y<3>>> break_points_y3
                = build_random_non_uniform_break_points(y3_min, y3_max, y3_ncells);

        std::vector<Coord<X<1>>> interpolation_points_x1 = break_points_x1;
        std::vector<Coord<X<2>>> interpolation_points_x2 = break_points_x2;
        std::vector<Coord<X<3>>> interpolation_points_x3 = break_points_x3;

        std::vector<Coord<Y<1>>> interpolation_points_y1 = break_points_y1;
        std::vector<Coord<Y<2>>> interpolation_points_y2 = break_points_y2;
        std::vector<Coord<Y<3>>> interpolation_points_y3 = break_points_y3;

        // Patch 1 ...............................................................................
        ddc::init_discrete_space<BSplinesX<1>>(break_points_x1);
        ddc::init_discrete_space<BSplinesY<1>>(break_points_y1);

        ddc::init_discrete_space<GridX<1>>(interpolation_points_x1);
        ddc::init_discrete_space<GridY<1>>(interpolation_points_y1);

        // Patch 2 ...............................................................................
        ddc::init_discrete_space<BSplinesX<2>>(break_points_x2);
        ddc::init_discrete_space<BSplinesY<2>>(break_points_y2);

        ddc::init_discrete_space<GridX<2>>(interpolation_points_x2);
        ddc::init_discrete_space<GridY<2>>(interpolation_points_y2);

        // Patch 3 ...............................................................................
        ddc::init_discrete_space<BSplinesX<3>>(break_points_x3);
        ddc::init_discrete_space<BSplinesY<3>>(break_points_y3);

        ddc::init_discrete_space<GridX<3>>(interpolation_points_x3);
        ddc::init_discrete_space<GridY<3>>(interpolation_points_y3);


        // Equivalent global domain ..............................................................
        std::vector<Coord<Xg>> break_points_xg;
        std::vector<Coord<Xg>> interpolation_points_xg;

        break_points_x1.pop_back();
        interpolation_points_x1.pop_back();
        fill_in(break_points_xg, break_points_x1);
        fill_in(interpolation_points_xg, interpolation_points_x1);

        break_points_x2.pop_back();
        interpolation_points_x2.pop_back();
        fill_in(break_points_xg, break_points_x2);
        fill_in(interpolation_points_xg, interpolation_points_x2);

        fill_in(break_points_xg, break_points_x3);
        fill_in(interpolation_points_xg, interpolation_points_x3);


        std::vector<Coord<Yg>> break_points_yg;
        std::vector<Coord<Yg>> interpolation_points_yg;

        fill_in(break_points_yg, break_points_y1);
        fill_in(interpolation_points_yg, interpolation_points_y1);


        ddc::init_discrete_space<BSplinesXg>(break_points_xg);
        ddc::init_discrete_space<BSplinesYg>(break_points_yg);

        ddc::init_discrete_space<GridXg>(interpolation_points_xg);
        ddc::init_discrete_space<GridYg>(interpolation_points_yg);
    }



    // TEST OPERATORS ============================================================================
    /** @brief Check that the local grids and the equivalent global grid match together for a 
     * given patch. The integers x_shift and y_shift correspond to a shift in the indices to match
     * with the correct index on the global grid. 
     */
    template <class Patch>
    void check_interpolation_grids(
            typename Patch::IdxRange12 const& idx_range,
            int const x_shift,
            int const y_shift)
    {
        using Grid1 = typename Patch::Grid1;
        using Grid2 = typename Patch::Grid2;
        ddc::for_each(idx_range, [&](typename Patch::Idx12 const& idx) {
            IdxStep<Grid1> idx_x(Idx<Grid1>(idx) - IdxRange<Grid1>(idx_range).front());
            IdxStep<Grid2> idx_y(Idx<Grid2>(idx) - IdxRange<Grid2>(idx_range).front());
            Idx<GridXg, GridYg> idx_g(idx_x.value() + x_shift, idx_y.value() + y_shift);
            EXPECT_NEAR(
                    ddc::coordinate(Idx<Grid1>(idx)),
                    ddc::coordinate(Idx<GridXg>(idx_g)),
                    1e-15);
            EXPECT_NEAR(
                    ddc::coordinate(Idx<Grid2>(idx)),
                    ddc::coordinate(Idx<GridYg>(idx_g)),
                    1e-15);
        });
    }


    /** @brief Check agreement between the computed x-derivatives and the global x-derivatives at 
     * the interfaces. 
     */
    template <class... Patches>
    void check_all_x_derivatives(
            MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_range_slices_dx)
    {
        (check_x_derivatives<Patches>(
                 functions_and_derivs.template get<Patches>(),
                 evaluator_g,
                 function_g_coef,
                 typename Patches::IdxRange1(idx_ranges.template get<Patches>()),
                 idx_range_slices_dx.template get<Patches>()),
         ...);
    }

    /** @brief Check agreement between the computed x-derivatives and the global x-derivatives at 
     * the interfaces for a given patch.
     */
    template <class Patch>
    void check_x_derivatives(
            DerivFieldOnPatch_host<Patch> function_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            typename Patch::IdxRange1 const& idx_range_perp,
            IdxRange1SliceOnPatch<Patch> const& idx_range_slice_dx)
    {
        using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
        Idx<DerivX> idx_deriv(1);

        Idx<DerivX, typename Patch::Grid1> idx_slice_min(idx_deriv, idx_range_slice_dx.front());
        Idx<DerivX, typename Patch::Grid1> idx_slice_max(idx_deriv, idx_range_slice_dx.back());

        DField<IdxRange<typename Patch::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_xmin_extracted = function_and_derivs[idx_slice_min];
        DField<IdxRange<typename Patch::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_xmax_extracted = function_and_derivs[idx_slice_max];

        typename Patch::IdxRange2 idx_range_par(get_idx_range(derivs_xmin_extracted));
        ddc::for_each(idx_range_par, [&](typename Patch::Idx2 const& idx_par) {
            typename Patch::Idx12 idx_min(idx_range_perp.front(), idx_par);
            typename Patch::Idx12 idx_max(idx_range_perp.back(), idx_par);
            Coord<Xg, Yg> interface_coord_min(get_global_coord<Xg, Yg>(ddc::coordinate(idx_min)));
            Coord<Xg, Yg> interface_coord_max(get_global_coord<Xg, Yg>(ddc::coordinate(idx_max)));

            double const global_deriv_min
                    = evaluator_g.deriv_dim_1(interface_coord_min, function_g_coef);
            double const global_deriv_max
                    = evaluator_g.deriv_dim_1(interface_coord_max, function_g_coef);

            EXPECT_NEAR(derivs_xmin_extracted(idx_par), global_deriv_min, 5e-14);
            EXPECT_NEAR(derivs_xmax_extracted(idx_par), global_deriv_max, 5e-14);
        });
    }


    /** @brief Check agreement between the computed y-derivatives and the global y-derivatives at 
     * the interfaces.
     */
    template <class... Patches>
    void check_all_y_derivatives(
            MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_range_slices_dy)
    {
        (check_y_derivatives<Patches>(
                 functions_and_derivs.template get<Patches>(),
                 evaluator_g,
                 function_g_coef,
                 typename Patches::IdxRange2(idx_ranges.template get<Patches>()),
                 idx_range_slices_dy.template get<Patches>()),
         ...);
    }

    /** @brief Check agreement between the computed y-derivatives and the global y-derivatives at 
     * the interfaces for a given patch.
     */
    template <class Patch>
    void check_y_derivatives(
            DerivFieldOnPatch_host<Patch> function_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            typename Patch::IdxRange2 const& idx_range_perp,
            IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy)
    {
        using DerivY = typename ddc::Deriv<typename Patch::Dim2>;
        Idx<DerivY> idx_deriv(1);

        Idx<DerivY, typename Patch::Grid2> idx_slice_min(idx_deriv, idx_range_slice_dy.front());
        Idx<DerivY, typename Patch::Grid2> idx_slice_max(idx_deriv, idx_range_slice_dy.back());

        DField<IdxRange<typename Patch::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_ymin_extracted = function_and_derivs[idx_slice_min];
        DField<IdxRange<typename Patch::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
                derivs_ymax_extracted = function_and_derivs[idx_slice_max];

        typename Patch::IdxRange1 idx_range_par(get_idx_range(derivs_ymin_extracted));
        ddc::for_each(idx_range_par, [&](typename Patch::Idx1 const& idx_par) {
            typename Patch::Idx12 idx_min(idx_par, idx_range_perp.front());
            typename Patch::Idx12 idx_max(idx_par, idx_range_perp.back());
            Coord<Xg, Yg> interface_coord_min(get_global_coord<Xg, Yg>(ddc::coordinate(idx_min)));
            Coord<Xg, Yg> interface_coord_max(get_global_coord<Xg, Yg>(ddc::coordinate(idx_max)));

            double const global_deriv_min
                    = evaluator_g.deriv_dim_2(interface_coord_min, function_g_coef);
            double const global_deriv_max
                    = evaluator_g.deriv_dim_2(interface_coord_max, function_g_coef);

            EXPECT_NEAR(derivs_ymin_extracted(idx_par), global_deriv_min, 5e-14);
            EXPECT_NEAR(derivs_ymax_extracted(idx_par), global_deriv_max, 5e-14);
        });
    }


    /** @brief Check agreement between the computed cross-derivatives and the global cross-derivatives at 
     * the interfaces.
     */
    template <class... Patches>
    void check_all_xy_derivatives(
            MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_range_slices_dx,
            MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_range_slices_dy)
    {
        (check_xy_derivatives<Patches>(
                 functions_and_derivs.template get<Patches>(),
                 evaluator_g,
                 function_g_coef,
                 idx_ranges.template get<Patches>(),
                 idx_range_slices_dx.template get<Patches>(),
                 idx_range_slices_dy.template get<Patches>()),
         ...);
    }

    /** @brief Check agreement between the computed cross-derivatives and the global cross-derivatives at 
     * the interfaces for a given patch.
     */
    template <class Patch>
    void check_xy_derivatives(
            DerivFieldOnPatch_host<Patch> function_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
            typename Patch::IdxRange12 const& idx_range,
            IdxRange1SliceOnPatch<Patch> const& idx_range_slice_dx,
            IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy)
    {
        using GridX = typename Patch::Grid1;
        using GridY = typename Patch::Grid2;
        using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
        using DerivY = typename ddc::Deriv<typename Patch::Dim2>;

        Idx<DerivX> idx_deriv_x(1);
        Idx<DerivY> idx_deriv_y(1);

        Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_min_min(
                idx_deriv_x,
                idx_range_slice_dx.front(),
                idx_deriv_y,
                idx_range_slice_dy.front());

        Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_max_min(
                idx_deriv_x,
                idx_range_slice_dx.back(),
                idx_deriv_y,
                idx_range_slice_dy.front());

        Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_min_max(
                idx_deriv_x,
                idx_range_slice_dx.front(),
                idx_deriv_y,
                idx_range_slice_dy.back());

        Idx<DerivX, GridX, DerivY, GridY> idx_cross_deriv_max_max(
                idx_deriv_x,
                idx_range_slice_dx.back(),
                idx_deriv_y,
                idx_range_slice_dy.back());

        IdxRange<GridX> idx_range_x(idx_range);
        IdxRange<GridY> idx_range_y(idx_range);

        Idx<GridX, GridY> idx_min_min(idx_range_x.front(), idx_range_y.front());
        Idx<GridX, GridY> idx_max_min(idx_range_x.back(), idx_range_y.front());
        Idx<GridX, GridY> idx_min_max(idx_range_x.front(), idx_range_y.back());
        Idx<GridX, GridY> idx_max_max(idx_range_x.back(), idx_range_y.back());

        Coord<Xg, Yg> interface_coord_min_min(
                get_global_coord<Xg, Yg>(ddc::coordinate(idx_min_min)));
        Coord<Xg, Yg> interface_coord_max_min(
                get_global_coord<Xg, Yg>(ddc::coordinate(idx_max_min)));
        Coord<Xg, Yg> interface_coord_min_max(
                get_global_coord<Xg, Yg>(ddc::coordinate(idx_min_max)));
        Coord<Xg, Yg> interface_coord_max_max(
                get_global_coord<Xg, Yg>(ddc::coordinate(idx_max_max)));

        double const global_deriv_min_min
                = evaluator_g.deriv_1_and_2(interface_coord_min_min, function_g_coef);
        double const global_deriv_max_min
                = evaluator_g.deriv_1_and_2(interface_coord_max_min, function_g_coef);
        double const global_deriv_min_max
                = evaluator_g.deriv_1_and_2(interface_coord_min_max, function_g_coef);
        double const global_deriv_max_max
                = evaluator_g.deriv_1_and_2(interface_coord_max_max, function_g_coef);

        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_min_min), global_deriv_min_min, 5e-13);
        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_max_min), global_deriv_max_min, 5e-13);
        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_min_max), global_deriv_min_max, 5e-13);
        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_max_max), global_deriv_max_max, 5e-13);
    }

    /** @brief Check agreement between the local splines defined with the computed derivatives 
     * and the global spline. 
     */
    template <class... Patches>
    void check_all_spline_representation_agreement(
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_range_slices_x,
            MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_range_slices_y,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
    {
        (check_spline_representation_agreement<Patches>(
                 idx_ranges.template get<Patches>(),
                 idx_range_slices_x.template get<Patches>(),
                 idx_range_slices_y.template get<Patches>(),
                 functions_and_derivs.template get<Patches>(),
                 evaluator_g,
                 function_g_coef),
         ...);
    };

    /** @brief Check agreement between the local splines defined with the computed derivatives 
     * and the global spline for a given patch.  
     */
    template <class Patch>
    void check_spline_representation_agreement(
            typename Patch::IdxRange12 const& idx_range_xy,
            IdxRange1SliceOnPatch<Patch> const& idx_range_slice_x,
            IdxRange2SliceOnPatch<Patch> const& idx_range_slice_y,
            DerivFieldOnPatch_host<Patch> function_and_derivs,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
    {
        // check_local_and_global_splines_agreement<
        //         Patch,
        //         ddc::BoundCond::HERMITE,
        //         ddc::BoundCond::HERMITE,
        //         ddc::BoundCond::HERMITE,
        //         ddc::BoundCond::HERMITE,
        //         BSplinesXg,
        //         BSplinesYg>(
        //         idx_range_xy,
        //         idx_range_slice_x,
        //         idx_range_slice_y,
        //         function_and_derivs,
        //         evaluator_g,
        //         function_g_coef);

        using DimX = typename Patch::Dim1;
        using DimY = typename Patch::Dim2;
        using DerivX = typename ddc::Deriv<DimX>;
        using DerivY = typename ddc::Deriv<DimY>;
        using GridX = typename Patch::Grid1;
        using GridY = typename Patch::Grid2;

        using Xg = typename BSplinesXg::continuous_dimension_type;
        using Yg = typename BSplinesYg::continuous_dimension_type;

        // Build the local spline representation -------------------------------------------------
        ddc::SplineBuilder2D<
                HostExecSpace,
                typename HostExecSpace::memory_space,
                typename Patch::BSplines1,
                typename Patch::BSplines2,
                typename Patch::Grid1,
                typename Patch::Grid2,
                ddc::BoundCond::HERMITE,
                ddc::BoundCond::HERMITE,
                ddc::BoundCond::HERMITE,
                ddc::BoundCond::HERMITE,
                ddc::SplineSolver::LAPACK>
                builder(idx_range_xy);

        SplineCoeffMemOnPatch_2D_host<Patch> function_coef_alloc(
                builder.batched_spline_domain(idx_range_xy));
        SplineCoeffOnPatch_2D_host<Patch> function_coef = get_field(function_coef_alloc);

        // Get the fields on the right layout.
        // --- extract each fields.
        // ------ function.
        DField<IdxRange<GridX, GridY>, Kokkos::HostSpace, Kokkos::layout_stride> function_extracted
                = function_and_derivs.get_values_field();

        // ------ first derivatives.
        IdxRange<DerivX> idx_range_deriv_x(Idx<DerivX>(1), IdxStep<DerivX>(1));
        IdxRange<DerivY> idx_range_deriv_y(Idx<DerivY>(1), IdxStep<DerivY>(1));
        IdxRange<DerivX, DerivY> idx_range_deriv_x_deriv_y(idx_range_deriv_x, idx_range_deriv_y);

        IdxRange<GridX> idx_range_x(idx_range_xy);
        IdxRange<GridY> idx_range_y(idx_range_xy);

        Idx<GridX> idx_slice_xmin(idx_range_slice_x.front());
        Idx<GridX> idx_slice_xmax(idx_range_slice_x.back());
        Idx<GridY> idx_slice_ymin(idx_range_slice_y.front());
        Idx<GridY> idx_slice_ymax(idx_range_slice_y.back());

        IdxRange<GridX> idx_range_slice_xmin(idx_slice_xmin, IdxStep<GridX>(1));
        IdxRange<GridX> idx_range_slice_xmax(idx_slice_xmax, IdxStep<GridX>(1));
        IdxRange<GridY> idx_range_slice_ymin(idx_slice_ymin, IdxStep<GridY>(1));
        IdxRange<GridY> idx_range_slice_ymax(idx_slice_ymax, IdxStep<GridY>(1));

        IdxRange<DerivX, GridX, GridY>
                idx_range_deriv_xmin(idx_range_deriv_x, idx_range_slice_xmin, idx_range_y);
        IdxRange<DerivX, GridX, GridY>
                idx_range_deriv_xmax(idx_range_deriv_x, idx_range_slice_xmax, idx_range_y);
        IdxRange<GridX, DerivY, GridY>
                idx_range_deriv_ymin(idx_range_x, idx_range_deriv_y, idx_range_slice_ymin);
        IdxRange<GridX, DerivY, GridY>
                idx_range_deriv_ymax(idx_range_x, idx_range_deriv_y, idx_range_slice_ymax);

        DField<IdxRange<DerivX, GridY>, Kokkos::HostSpace, Kokkos::layout_stride>
                deriv_xmin_extracted = function_and_derivs[idx_range_deriv_xmin][idx_slice_xmin];
        DField<IdxRange<DerivX, GridY>, Kokkos::HostSpace, Kokkos::layout_stride>
                deriv_xmax_extracted = function_and_derivs[idx_range_deriv_xmax][idx_slice_xmax];

        DField<IdxRange<GridX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
                deriv_ymin_extracted = function_and_derivs[idx_range_deriv_ymin][idx_slice_ymin];
        DField<IdxRange<GridX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
                deriv_ymax_extracted = function_and_derivs[idx_range_deriv_ymax][idx_slice_ymax];

        // ------ cross-derivatives.
        IdxRange<DerivX, GridX, DerivY, GridY> idx_range_deriv_xy_min_min(
                idx_range_deriv_x,
                idx_range_slice_xmin,
                idx_range_deriv_y,
                idx_range_slice_ymin);
        IdxRange<DerivX, GridX, DerivY, GridY> idx_range_deriv_xy_max_min(
                idx_range_deriv_x,
                idx_range_slice_xmax,
                idx_range_deriv_y,
                idx_range_slice_ymin);
        IdxRange<DerivX, GridX, DerivY, GridY> idx_range_deriv_xy_min_max(
                idx_range_deriv_x,
                idx_range_slice_xmin,
                idx_range_deriv_y,
                idx_range_slice_ymax);
        IdxRange<DerivX, GridX, DerivY, GridY> idx_range_deriv_xy_max_max(
                idx_range_deriv_x,
                idx_range_slice_xmax,
                idx_range_deriv_y,
                idx_range_slice_ymax);

        DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
                deriv_xy_min_min_extracted
                = function_and_derivs[idx_range_deriv_xy_min_min][idx_slice_xmin][idx_slice_ymin];
        DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
                deriv_xy_max_min_extracted
                = function_and_derivs[idx_range_deriv_xy_max_min][idx_slice_xmax][idx_slice_ymin];
        DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
                deriv_xy_min_max_extracted
                = function_and_derivs[idx_range_deriv_xy_min_max][idx_slice_xmin][idx_slice_ymax];
        DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
                deriv_xy_max_max_extracted
                = function_and_derivs[idx_range_deriv_xy_max_max][idx_slice_xmax][idx_slice_ymax];


        // --- define fields on the correct layout.
        // ------ allocate memory
        host_t<DFieldMem<IdxRange<GridX, GridY>>> function_alloc(idx_range_xy);

        IdxRange<DerivX, GridY> idx_range_deriv_x_y(idx_range_deriv_x, idx_range_y);
        host_t<DFieldMem<IdxRange<DerivX, GridY>>> deriv_xmin_alloc(idx_range_deriv_x_y);
        host_t<DFieldMem<IdxRange<DerivX, GridY>>> deriv_xmax_alloc(idx_range_deriv_x_y);

        IdxRange<GridX, DerivY> idx_range_x_deriv_y(idx_range_x, idx_range_deriv_y);
        host_t<DFieldMem<IdxRange<GridX, DerivY>>> deriv_ymin_alloc(idx_range_x_deriv_y);
        host_t<DFieldMem<IdxRange<GridX, DerivY>>> deriv_ymax_alloc(idx_range_x_deriv_y);

        host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_min_min_alloc(
                idx_range_deriv_x_deriv_y);
        host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_max_min_alloc(
                idx_range_deriv_x_deriv_y);
        host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_min_max_alloc(
                idx_range_deriv_x_deriv_y);
        host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_max_max_alloc(
                idx_range_deriv_x_deriv_y);

        // ------ define span
        host_t<DField<IdxRange<GridX, GridY>>> function(function_alloc);

        host_t<DField<IdxRange<DerivX, GridY>>> deriv_xmin(deriv_xmin_alloc);
        host_t<DField<IdxRange<DerivX, GridY>>> deriv_xmax(deriv_xmax_alloc);

        host_t<DField<IdxRange<GridX, DerivY>>> deriv_ymin(deriv_ymin_alloc);
        host_t<DField<IdxRange<GridX, DerivY>>> deriv_ymax(deriv_ymax_alloc);

        host_t<DField<IdxRange<DerivX, DerivY>>> deriv_xy_min_min(deriv_xy_min_min_alloc);
        host_t<DField<IdxRange<DerivX, DerivY>>> deriv_xy_max_min(deriv_xy_max_min_alloc);
        host_t<DField<IdxRange<DerivX, DerivY>>> deriv_xy_min_max(deriv_xy_min_max_alloc);
        host_t<DField<IdxRange<DerivX, DerivY>>> deriv_xy_max_max(deriv_xy_max_max_alloc);

        // ------ initialise data from the fields on layout stride.
        ddc::for_each(idx_range_xy, [&](Idx<GridX, GridY> const idx) {
            function(idx) = function_extracted(idx);
        });

        ddc::for_each(idx_range_deriv_x_y, [&](Idx<DerivX, GridY> const idx) {
            deriv_xmin(idx) = deriv_xmin_extracted(idx);
            deriv_xmax(idx) = deriv_xmax_extracted(idx);
        });

        ddc::for_each(idx_range_x_deriv_y, [&](Idx<GridX, DerivY> const idx) {
            deriv_ymin(idx) = deriv_ymin_extracted(idx);
            deriv_ymax(idx) = deriv_ymax_extracted(idx);
        });

        Idx<DerivX, DerivY> const idx_cross_deriv(idx_range_deriv_x_deriv_y.front());
        deriv_xy_min_min(idx_cross_deriv) = deriv_xy_min_min_extracted(idx_cross_deriv);
        deriv_xy_max_min(idx_cross_deriv) = deriv_xy_max_min_extracted(idx_cross_deriv);
        deriv_xy_min_max(idx_cross_deriv) = deriv_xy_min_max_extracted(idx_cross_deriv);
        deriv_xy_max_max(idx_cross_deriv) = deriv_xy_max_max_extracted(idx_cross_deriv);


        builder(function_coef,
                get_const_field(function),
                std::optional(get_const_field(deriv_xmin)),
                std::optional(get_const_field(deriv_xmax)),
                std::optional(get_const_field(deriv_ymin)),
                std::optional(get_const_field(deriv_ymax)),
                std::optional(get_const_field(deriv_xy_min_min)),
                std::optional(get_const_field(deriv_xy_max_min)),
                std::optional(get_const_field(deriv_xy_min_max)),
                std::optional(get_const_field(deriv_xy_max_max)));

        // Define local spline evaluator ---------------------------------------------------------
        Coord<DimX> const x_min(ddc::discrete_space<typename Patch::BSplines1>().rmin());
        Coord<DimX> const x_max(ddc::discrete_space<typename Patch::BSplines1>().rmax());
        Coord<DimY> const y_min(ddc::discrete_space<typename Patch::BSplines2>().rmin());
        Coord<DimY> const y_max(ddc::discrete_space<typename Patch::BSplines2>().rmax());
        ddc::ConstantExtrapolationRule<DimX, DimY> bc_xmin(x_min, y_min, y_max);
        ddc::ConstantExtrapolationRule<DimX, DimY> bc_xmax(x_max, y_min, y_max);
        ddc::ConstantExtrapolationRule<DimY, DimX> bc_ymin(y_min, x_min, x_max);
        ddc::ConstantExtrapolationRule<DimY, DimX> bc_ymax(y_max, x_min, x_max);
        ddc::SplineEvaluator2D<
                HostExecSpace,
                typename HostExecSpace::memory_space,
                typename Patch::BSplines1,
                typename Patch::BSplines2,
                typename Patch::Grid1,
                typename Patch::Grid2,
                ddc::ConstantExtrapolationRule<DimX, DimY>,
                ddc::ConstantExtrapolationRule<DimX, DimY>,
                ddc::ConstantExtrapolationRule<DimY, DimX>,
                ddc::ConstantExtrapolationRule<DimY, DimX>>
                evaluator(bc_xmin, bc_xmax, bc_ymin, bc_ymax);

        // Define evaluation points at the centre of the cells -----------------------------------
        host_t<CoordFieldMemOnPatch<Patch>> eval_points_alloc(idx_range_xy);
        host_t<CoordFieldOnPatch<Patch>> eval_points(eval_points_alloc);

        // The evaluation points are placed in the middle of the cells, except for the last one.
        ddc::for_each(idx_range_xy, [&](typename Patch::Idx12 const idx) {
            Coord<DimX, DimY> const mesh_point(ddc::coordinate(idx));
            typename Patch::Idx1 idx_1(idx);
            typename Patch::Idx2 idx_2(idx);
            Coord<DimX> dx = (idx_1 != typename Patch::IdxRange1(idx_range_xy).back())
                             * distance_at_right(idx_1);
            Coord<DimY> dy = (idx_2 != typename Patch::IdxRange2(idx_range_xy).back())
                             * distance_at_right(idx_2);
            eval_points(idx) = mesh_point + Coord<DimX, DimY>(dx, dy);
        });

        // Evaluate and compare the local and global spline representations ----------------------
        ddc::for_each(idx_range_xy, [&](typename Patch::Idx12 const idx) {
            Coord<Xg, Yg> const eval_point_g(get_global_coord<Xg, Yg>(eval_points(idx)));

            double local_spline = evaluator(eval_points(idx), get_const_field(function_coef));
            double global_spline = evaluator_g(eval_point_g, get_const_field(function_g_coef));

            EXPECT_NEAR(local_spline, global_spline, 1e-14);
        });
    }
};

} // end namespace



// // Check that the local grids and the equivalent global grid match together.
TEST_F(InterfaceDerivativeMatrixHermiteTest, InterpolationPointsCheck)
{
    int const x_shift1 = x1_ncells.value();
    int const x_shift2 = x1_ncells.value() + x2_ncells.value();

    check_interpolation_grids<Patch1>(idx_range_xy1, 0, 0);
    check_interpolation_grids<Patch2>(idx_range_xy2, x_shift1, 0);
    check_interpolation_grids<Patch3>(idx_range_xy3, x_shift2, 0);
}



/*
 * Here, we consider only the 1 | 2 | 3 strip. We assume that the equivalant global 
 * domain is not periodic but the derivatives at the boundaries are known. 
 * It allows us to test the InterfaceDerivativeMatrix for ddc::BoundCond::HERMITE.
 */
TEST_F(InterfaceDerivativeMatrixHermiteTest, CheckForHermiteBc)
{
    // Instantiate the derivatives calculators ---------------------------------------------------
    // SingleInterfaceDerivativesCalculators for interfaces along y (periodic).
    SingleInterfaceDerivativesCalculator<Interface_1_2> const
            derivatives_calculator_1_2(idx_range_xy1, idx_range_xy2);
    SingleInterfaceDerivativesCalculator<Interface_2_3> const
            derivatives_calculator_2_3(idx_range_xy2, idx_range_xy3);

    // Collect the derivative calculators --------------------------------------------------------
    SingleInterfaceDerivativesCalculatorCollection
            deriv_calculators_collect(derivatives_calculator_1_2, derivatives_calculator_2_3);

    // Collect the index ranges ------------------------------------------------------------------
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2, Patch3>
            idx_ranges(idx_range_xy1, idx_range_xy2, idx_range_xy3);

    // Instantiate the matrix calculators --------------------------------------------------------
    InterfaceDerivativeMatrix<
            Connectivity,
            GridX<1>,
            decltype(deriv_calculators_collect),
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::HERMITE,
            Kokkos::DefaultHostExecutionSpace,
            Patch1,
            Patch2,
            Patch3>
            matrix(idx_ranges, deriv_calculators_collect);


    // Instantiate DerivField ====================================================================
    // Instantiate index range slices ------------------------------------------------------------
    IdxRangeSlice<GridX<1>>
            idx_range_slice_dx1(idx_range_x1.front(), IdxStep<GridX<1>>(2), idx_range_x1.extents());
    IdxRangeSlice<GridX<2>>
            idx_range_slice_dx2(idx_range_x2.front(), IdxStep<GridX<2>>(2), idx_range_x2.extents());
    IdxRangeSlice<GridX<3>>
            idx_range_slice_dx3(idx_range_x3.front(), IdxStep<GridX<3>>(2), idx_range_x3.extents());

    IdxRangeSlice<GridY<1>>
            idx_range_slice_dy1(idx_range_y1.front(), IdxStep<GridY<1>>(2), idx_range_y1.extents());
    IdxRangeSlice<GridY<2>>
            idx_range_slice_dy2(idx_range_y2.front(), IdxStep<GridY<2>>(2), idx_range_y2.extents());
    IdxRangeSlice<GridY<3>>
            idx_range_slice_dy3(idx_range_y3.front(), IdxStep<GridY<3>>(2), idx_range_y3.extents());

    // Collect the index range slices.
    MultipatchType<IdxRange1SliceOnPatch, Patch1, Patch2, Patch3>
            idx_ranges_slice_dx(idx_range_slice_dx1, idx_range_slice_dx2, idx_range_slice_dx3);

    MultipatchType<IdxRange2SliceOnPatch, Patch1, Patch2, Patch3>
            idx_ranges_slice_dy(idx_range_slice_dy1, idx_range_slice_dy2, idx_range_slice_dy3);

    // Instantiate DerivField --------------------------------------------------------------------
    DerivFieldMemOnPatch_host<Patch1>
            function_and_derivs_1_alloc(idx_range_xy1, idx_range_slice_dx1, idx_range_slice_dy1);
    DerivFieldMemOnPatch_host<Patch2>
            function_and_derivs_2_alloc(idx_range_xy2, idx_range_slice_dx2, idx_range_slice_dy2);
    DerivFieldMemOnPatch_host<Patch3>
            function_and_derivs_3_alloc(idx_range_xy3, idx_range_slice_dx3, idx_range_slice_dy3);

    DerivFieldOnPatch_host<Patch1> function_and_derivs_1(function_and_derivs_1_alloc);
    DerivFieldOnPatch_host<Patch2> function_and_derivs_2(function_and_derivs_2_alloc);
    DerivFieldOnPatch_host<Patch3> function_and_derivs_3(function_and_derivs_3_alloc);

    // Collect the fields with derivatives.
    MultipatchField<DerivFieldOnPatch_host, Patch1, Patch2, Patch3> functions_and_derivs(
            function_and_derivs_1,
            function_and_derivs_2,
            function_and_derivs_3);

    // Instantiate the field of global function values. No derivatives needed.
    host_t<DFieldMem<IdxRange<GridXg, GridYg>>> function_g_alloc(idx_range_xy_g);
    host_t<DField<IdxRange<GridXg, GridYg>>> function_g = get_field(function_g_alloc);

    // Initialise the data =======================================================================
    // --- the function values.
    initialise_all_functions(functions_and_derivs);
    initialise_2D_function<GridXg, GridYg>(function_g);

    // --- the derivatives
    Idx<DerivXg> first_dxg(1);
    IdxStep<DerivXg> n_deriv_xg(1);
    IdxRange<DerivXg> idx_range_deriv_xg(first_dxg, n_deriv_xg);

    Idx<DerivYg> first_dyg(1);
    IdxStep<DerivYg> n_deriv_yg(1);
    IdxRange<DerivYg> idx_range_deriv_yg(first_dyg, n_deriv_yg);

    IdxRange<DerivXg, GridYg> idx_range_dxg_yg(idx_range_deriv_xg, idx_range_yg);
    IdxRange<GridXg, DerivYg> idx_range_xg_dyg(idx_range_xg, idx_range_deriv_yg);
    IdxRange<DerivXg, DerivYg> idx_range_dxg_dyg(idx_range_deriv_xg, idx_range_deriv_yg);

    host_t<DFieldMem<IdxRange<DerivXg, GridYg>>> derivs_xgmin_alloc(idx_range_dxg_yg);
    host_t<DFieldMem<IdxRange<DerivXg, GridYg>>> derivs_xgmax_alloc(idx_range_dxg_yg);
    host_t<DFieldMem<IdxRange<GridXg, DerivYg>>> derivs_ygmin_alloc(idx_range_xg_dyg);
    host_t<DFieldMem<IdxRange<GridXg, DerivYg>>> derivs_ygmax_alloc(idx_range_xg_dyg);

    host_t<DFieldMem<IdxRange<DerivXg, DerivYg>>> derivs_xyg_min_min_alloc(idx_range_dxg_dyg);
    host_t<DFieldMem<IdxRange<DerivXg, DerivYg>>> derivs_xyg_max_min_alloc(idx_range_dxg_dyg);
    host_t<DFieldMem<IdxRange<DerivXg, DerivYg>>> derivs_xyg_min_max_alloc(idx_range_dxg_dyg);
    host_t<DFieldMem<IdxRange<DerivXg, DerivYg>>> derivs_xyg_max_max_alloc(idx_range_dxg_dyg);

    host_t<DField<IdxRange<DerivXg, GridYg>>> derivs_xgmin = get_field(derivs_xgmin_alloc);
    host_t<DField<IdxRange<DerivXg, GridYg>>> derivs_xgmax = get_field(derivs_xgmax_alloc);
    host_t<DField<IdxRange<GridXg, DerivYg>>> derivs_ygmin = get_field(derivs_ygmin_alloc);
    host_t<DField<IdxRange<GridXg, DerivYg>>> derivs_ygmax = get_field(derivs_ygmax_alloc);

    host_t<DField<IdxRange<DerivXg, DerivYg>>> derivs_xyg_min_min(derivs_xyg_min_min_alloc);
    host_t<DField<IdxRange<DerivXg, DerivYg>>> derivs_xyg_max_min(derivs_xyg_max_min_alloc);
    host_t<DField<IdxRange<DerivXg, DerivYg>>> derivs_xyg_min_max(derivs_xyg_min_max_alloc);
    host_t<DField<IdxRange<DerivXg, DerivYg>>> derivs_xyg_max_max(derivs_xyg_max_max_alloc);

    ddc::for_each(idx_range_dxg_yg, [&](Idx<DerivXg, GridYg> const idx) {
        double const xgmin = xg_min;
        double const xgmax = xg_max;
        double const yg = ddc::coordinate(Idx<GridYg>(idx));
        derivs_xgmin(idx) = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xgmin) * std::sin(yg);
        derivs_xgmax(idx) = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xgmax) * std::sin(yg);
    });
    ddc::for_each(idx_range_xg_dyg, [&](Idx<GridXg, DerivYg> const idx) {
        double const ygmin = yg_min;
        double const ygmax = yg_max;
        double const xg = ddc::coordinate(Idx<GridXg>(idx));
        derivs_ygmin(idx) = std::cos(2. / 3 * M_PI * xg) * std ::cos(ygmin);
        derivs_ygmax(idx) = std::cos(2. / 3 * M_PI * xg) * std ::cos(ygmax);
    });
    derivs_xyg_min_min(first_dxg, first_dyg)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xg_min) * std::sin(yg_min);
    derivs_xyg_max_min(first_dxg, first_dyg)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xg_max) * std::sin(yg_min);
    derivs_xyg_min_max(first_dxg, first_dyg)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xg_min) * std::sin(yg_max);
    derivs_xyg_max_max(first_dxg, first_dyg)
            = -2. / 3 * M_PI * std::sin(2. / 3 * M_PI * xg_max) * std::sin(yg_max);

    // --- the derivatives from an equivalent global spline.
    // ------- build global spline representation
    SplineRThetagBuilder builder_g(idx_range_xy_g);

    host_t<DFieldMem<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef_alloc(
            builder_g.batched_spline_domain(idx_range_xy_g));
    host_t<DField<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef
            = get_field(function_g_coef_alloc);

    builder_g(
            function_g_coef,
            get_const_field(function_g),
            std::optional(get_const_field(derivs_xgmin)),
            std::optional(get_const_field(derivs_xgmax)),
            std::optional(get_const_field(derivs_ygmin)),
            std::optional(get_const_field(derivs_ygmax)),
            std::optional(get_const_field(derivs_xyg_min_min)),
            std::optional(get_const_field(derivs_xyg_max_min)),
            std::optional(get_const_field(derivs_xyg_min_max)),
            std::optional(get_const_field(derivs_xyg_max_max)));

    host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const_function_g_coef
            = get_const_field(function_g_coef);

    // ------ global spline evaluator
    ddc::ConstantExtrapolationRule<Yg, Xg> bc_ymin_g(yg_min, xg_min, xg_max);
    ddc::ConstantExtrapolationRule<Yg, Xg> bc_ymax_g(yg_max, xg_min, xg_max);
    ddc::ConstantExtrapolationRule<Xg, Yg> bc_xmin_g(xg_min, yg_min, yg_max);
    ddc::ConstantExtrapolationRule<Xg, Yg> bc_xmax_g(xg_max, yg_min, yg_max);
    SplineRThetagEvaluator evaluator_g(bc_xmin_g, bc_xmax_g, bc_ymin_g, bc_ymax_g);

    // ------ intialise the first derivatives
    Idx<ddc::Deriv<X<1>>, GridX<1>>
            idx_slice_xmin_1(Idx<ddc::Deriv<X<1>>>(1), idx_range_slice_dx1.front());
    Idx<ddc::Deriv<X<3>>, GridX<3>>
            idx_slice_xmax_3(Idx<ddc::Deriv<X<3>>>(1), idx_range_slice_dx3.back());

    DField<IdxRange<typename Patch1::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_xmin_extracted = function_and_derivs_1[idx_slice_xmin_1];
    ddc::for_each(idx_range_y1, [&](Idx<GridY<1>> const& idx_par) {
        Idx<GridX<1>, GridY<1>> idx_min(idx_range_x1.front(), idx_par);
        Coord<Xg, Yg> interface_coord_min(get_global_coord<Xg, Yg>(ddc::coordinate(idx_min)));
        derivs_xmin_extracted(idx_par)
                = evaluator_g.deriv_dim_1(interface_coord_min, const_function_g_coef);
    });

    DField<IdxRange<typename Patch3::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_xmax_extracted = function_and_derivs_3[idx_slice_xmax_3];
    ddc::for_each(idx_range_y3, [&](Idx<GridY<3>> const& idx_par) {
        Idx<GridX<3>, GridY<3>> idx_min(idx_range_x3.back(), idx_par);
        Coord<Xg, Yg> interface_coord_min(get_global_coord<Xg, Yg>(ddc::coordinate(idx_min)));
        derivs_xmax_extracted(idx_par)
                = evaluator_g.deriv_dim_1(interface_coord_min, const_function_g_coef);
    });

    Idx<ddc::Deriv<Y<1>>, GridY<1>>
            idx_slice_ymin_1(Idx<ddc::Deriv<Y<1>>>(1), idx_range_slice_dy1.front());
    Idx<ddc::Deriv<Y<1>>, GridY<1>>
            idx_slice_ymax_1(Idx<ddc::Deriv<Y<1>>>(1), idx_range_slice_dy1.back());

    Idx<ddc::Deriv<Y<2>>, GridY<2>>
            idx_slice_ymin_2(Idx<ddc::Deriv<Y<2>>>(1), idx_range_slice_dy2.front());
    Idx<ddc::Deriv<Y<2>>, GridY<2>>
            idx_slice_ymax_2(Idx<ddc::Deriv<Y<2>>>(1), idx_range_slice_dy2.back());

    Idx<ddc::Deriv<Y<3>>, GridY<3>>
            idx_slice_ymin_3(Idx<ddc::Deriv<Y<3>>>(1), idx_range_slice_dy3.front());
    Idx<ddc::Deriv<Y<3>>, GridY<3>>
            idx_slice_ymax_3(Idx<ddc::Deriv<Y<3>>>(1), idx_range_slice_dy3.back());


    DField<IdxRange<typename Patch1::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymin_extracted_1 = function_and_derivs_1[idx_slice_ymin_1];
    DField<IdxRange<typename Patch1::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymax_extracted_1 = function_and_derivs_1[idx_slice_ymax_1];
    ddc::for_each(idx_range_x1, [&](Idx<GridX<1>> const& idx_par) {
        Idx<GridX<1>, GridY<1>> idx_min(idx_par, idx_range_y1.front());
        Idx<GridX<1>, GridY<1>> idx_max(idx_par, idx_range_y1.back());
        Coord<Xg, Yg> interface_coord_min(get_global_coord<Xg, Yg>(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(get_global_coord<Xg, Yg>(ddc::coordinate(idx_max)));
        derivs_ymin_extracted_1(idx_par)
                = evaluator_g.deriv_dim_2(interface_coord_min, const_function_g_coef);
        derivs_ymax_extracted_1(idx_par)
                = evaluator_g.deriv_dim_2(interface_coord_max, const_function_g_coef);
    });

    DField<IdxRange<typename Patch2::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymin_extracted_2 = function_and_derivs_2[idx_slice_ymin_2];
    DField<IdxRange<typename Patch2::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymax_extracted_2 = function_and_derivs_2[idx_slice_ymax_2];
    ddc::for_each(idx_range_x2, [&](Idx<GridX<2>> const& idx_par) {
        Idx<GridX<2>, GridY<2>> idx_min(idx_par, idx_range_y2.front());
        Idx<GridX<2>, GridY<2>> idx_max(idx_par, idx_range_y2.back());
        Coord<Xg, Yg> interface_coord_min(get_global_coord<Xg, Yg>(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(get_global_coord<Xg, Yg>(ddc::coordinate(idx_max)));
        derivs_ymin_extracted_2(idx_par)
                = evaluator_g.deriv_dim_2(interface_coord_min, const_function_g_coef);
        derivs_ymax_extracted_2(idx_par)
                = evaluator_g.deriv_dim_2(interface_coord_max, const_function_g_coef);
    });

    DField<IdxRange<typename Patch3::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymin_extracted_3 = function_and_derivs_3[idx_slice_ymin_3];
    DField<IdxRange<typename Patch3::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymax_extracted_3 = function_and_derivs_3[idx_slice_ymax_3];
    ddc::for_each(idx_range_x3, [&](Idx<GridX<3>> const& idx_par) {
        Idx<GridX<3>, GridY<3>> idx_min(idx_par, idx_range_y3.front());
        Idx<GridX<3>, GridY<3>> idx_max(idx_par, idx_range_y3.back());
        Coord<Xg, Yg> interface_coord_min(get_global_coord<Xg, Yg>(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(get_global_coord<Xg, Yg>(ddc::coordinate(idx_max)));
        derivs_ymin_extracted_3(idx_par)
                = evaluator_g.deriv_dim_2(interface_coord_min, const_function_g_coef);
        derivs_ymax_extracted_3(idx_par)
                = evaluator_g.deriv_dim_2(interface_coord_max, const_function_g_coef);
    });

    // ------ intialise the cross-derivatives
    Idx<ddc::Deriv<X<1>>, GridX<1>, ddc::Deriv<Y<1>>, GridY<1>> idx_cross_deriv_min_min_1(
            Idx<ddc::Deriv<X<1>>>(1),
            idx_range_slice_dx1.front(),
            Idx<ddc::Deriv<Y<1>>>(1),
            idx_range_slice_dy1.front());
    Idx<ddc::Deriv<X<1>>, GridX<1>, ddc::Deriv<Y<1>>, GridY<1>> idx_cross_deriv_max_min_1(
            Idx<ddc::Deriv<X<1>>>(1),
            idx_range_slice_dx1.back(),
            Idx<ddc::Deriv<Y<1>>>(1),
            idx_range_slice_dy1.front());
    Idx<ddc::Deriv<X<1>>, GridX<1>, ddc::Deriv<Y<1>>, GridY<1>> idx_cross_deriv_min_max_1(
            Idx<ddc::Deriv<X<1>>>(1),
            idx_range_slice_dx1.front(),
            Idx<ddc::Deriv<Y<1>>>(1),
            idx_range_slice_dy1.back());
    Idx<ddc::Deriv<X<1>>, GridX<1>, ddc::Deriv<Y<1>>, GridY<1>> idx_cross_deriv_max_max_1(
            Idx<ddc::Deriv<X<1>>>(1),
            idx_range_slice_dx1.back(),
            Idx<ddc::Deriv<Y<1>>>(1),
            idx_range_slice_dy1.back());

    Idx<ddc::Deriv<X<2>>, GridX<2>, ddc::Deriv<Y<2>>, GridY<2>> idx_cross_deriv_min_min_2(
            Idx<ddc::Deriv<X<2>>>(1),
            idx_range_slice_dx2.front(),
            Idx<ddc::Deriv<Y<2>>>(1),
            idx_range_slice_dy2.front());
    Idx<ddc::Deriv<X<2>>, GridX<2>, ddc::Deriv<Y<2>>, GridY<2>> idx_cross_deriv_max_min_2(
            Idx<ddc::Deriv<X<2>>>(1),
            idx_range_slice_dx2.back(),
            Idx<ddc::Deriv<Y<2>>>(1),
            idx_range_slice_dy2.front());
    Idx<ddc::Deriv<X<2>>, GridX<2>, ddc::Deriv<Y<2>>, GridY<2>> idx_cross_deriv_min_max_2(
            Idx<ddc::Deriv<X<2>>>(1),
            idx_range_slice_dx2.front(),
            Idx<ddc::Deriv<Y<2>>>(1),
            idx_range_slice_dy2.back());
    Idx<ddc::Deriv<X<2>>, GridX<2>, ddc::Deriv<Y<2>>, GridY<2>> idx_cross_deriv_max_max_2(
            Idx<ddc::Deriv<X<2>>>(1),
            idx_range_slice_dx2.back(),
            Idx<ddc::Deriv<Y<2>>>(1),
            idx_range_slice_dy2.back());

    Idx<ddc::Deriv<X<3>>, GridX<3>, ddc::Deriv<Y<3>>, GridY<3>> idx_cross_deriv_min_min_3(
            Idx<ddc::Deriv<X<3>>>(1),
            idx_range_slice_dx3.front(),
            Idx<ddc::Deriv<Y<3>>>(1),
            idx_range_slice_dy3.front());
    Idx<ddc::Deriv<X<3>>, GridX<3>, ddc::Deriv<Y<3>>, GridY<3>> idx_cross_deriv_max_min_3(
            Idx<ddc::Deriv<X<3>>>(1),
            idx_range_slice_dx3.back(),
            Idx<ddc::Deriv<Y<3>>>(1),
            idx_range_slice_dy3.front());
    Idx<ddc::Deriv<X<3>>, GridX<3>, ddc::Deriv<Y<3>>, GridY<3>> idx_cross_deriv_min_max_3(
            Idx<ddc::Deriv<X<3>>>(1),
            idx_range_slice_dx3.front(),
            Idx<ddc::Deriv<Y<3>>>(1),
            idx_range_slice_dy3.back());
    Idx<ddc::Deriv<X<3>>, GridX<3>, ddc::Deriv<Y<3>>, GridY<3>> idx_cross_deriv_max_max_3(
            Idx<ddc::Deriv<X<3>>>(1),
            idx_range_slice_dx3.back(),
            Idx<ddc::Deriv<Y<3>>>(1),
            idx_range_slice_dy3.back());

    function_and_derivs_1(idx_cross_deriv_min_min_1)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(0, 0), const_function_g_coef);
    function_and_derivs_1(idx_cross_deriv_max_min_1)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(1, 0), const_function_g_coef);
    function_and_derivs_1(idx_cross_deriv_min_max_1)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(0, 1), const_function_g_coef);
    function_and_derivs_1(idx_cross_deriv_max_max_1)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(1, 1), const_function_g_coef);

    function_and_derivs_2(idx_cross_deriv_min_min_2)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(1, 0), const_function_g_coef);
    function_and_derivs_2(idx_cross_deriv_max_min_2)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(2, 0), const_function_g_coef);
    function_and_derivs_2(idx_cross_deriv_min_max_2)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(1, 1), const_function_g_coef);
    function_and_derivs_2(idx_cross_deriv_max_max_2)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(2, 1), const_function_g_coef);

    function_and_derivs_3(idx_cross_deriv_min_min_3)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(2, 0), const_function_g_coef);
    function_and_derivs_3(idx_cross_deriv_max_min_3)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(3, 0), const_function_g_coef);
    function_and_derivs_3(idx_cross_deriv_min_max_3)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(2, 1), const_function_g_coef);
    function_and_derivs_3(idx_cross_deriv_max_max_3)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(3, 1), const_function_g_coef);


    // --- the first derivatives from the function values.
    matrix.solve_deriv(functions_and_derivs);

    // --- the cross-derivatives from the first derivatives.
    matrix.solve_cross_deriv(functions_and_derivs);

    // Test the values of the derivatives ========================================================
    // Check each derivatives ---
    check_all_x_derivatives(
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            idx_ranges,
            idx_ranges_slice_dx);

    check_all_y_derivatives(
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            idx_ranges,
            idx_ranges_slice_dy);

    check_all_xy_derivatives(
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            idx_ranges,
            idx_ranges_slice_dx,
            idx_ranges_slice_dy);

    // Check the whole spline representations ---
    check_all_spline_representation_agreement(
            idx_ranges,
            idx_ranges_slice_dx,
            idx_ranges_slice_dy,
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef);


    std::cout << "End of the test." << std::endl;
}
