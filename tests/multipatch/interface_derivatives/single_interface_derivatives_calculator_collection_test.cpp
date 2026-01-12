// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "3patches_2d_non_periodic_non_uniform.hpp"
#include "interface.hpp"
#include "interface_derivatives_test_utils.hpp"
#include "mesh_builder.hpp"
#include "multipatch_field.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "single_interface_derivatives_calculator.hpp"
#include "single_interface_derivatives_calculator_collection.hpp"
#include "view.hpp"


/*
    Test SingleInterfaceDerivativesCalculatorCollection on the following geometry:

        |  1  |  2  |  3  |

        with the global X and Y dimensions with additional points as closure condition 
        (ddc::BoundCond::GREVILLE).

*/

namespace {
// Multi-patch tags ---
using namespace non_periodic_non_uniform_2d_3patches;

// INTERFACES ------------------------------------------------------------------------------------
using Interface_12 = Interface<EastEdge<1>, WestEdge<2>, true>;
using Interface_23 = Interface<EastEdge<2>, WestEdge<3>, true>;


using HostExecSpace = Kokkos::DefaultHostExecutionSpace;


// Interpolation points type for the patches.
template <std::size_t PatchIdx, ddc::BoundCond BC1, ddc::BoundCond BC2>
using SplineInterpPointsX = ddcHelper::NonUniformInterpolationPoints<BSplinesX<PatchIdx>, BC1, BC2>;

template <std::size_t PatchIdx, ddc::BoundCond BC1, ddc::BoundCond BC2>
using SplineInterpPointsY = ddcHelper::NonUniformInterpolationPoints<BSplinesY<PatchIdx>, BC1, BC2>;


static constexpr ddc::BoundCond BCH = ddc::BoundCond::HERMITE;
static constexpr ddc::BoundCond BCG = ddc::BoundCond::GREVILLE;


// USEFUL FUNCTIONS ==============================================================================
/// @brief Initialise the function with f(x,y) = cos(2/3*pi*x)sin(y -0.3).
template <class Grid1, class Grid2>
void initialise_2D_function(host_t<DField<IdxRange<Grid1, Grid2>>> function)
{
    ddc::host_for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        double const xg = ddc::coordinate(Idx<Grid1>(idx));
        double const yg = ddc::coordinate(Idx<Grid2>(idx));
        function(idx) = Kokkos::cos(xg * 2. / 3. * M_PI) * Kokkos::sin(yg - 0.3);
    });
}

/// @brief Initialise all the functions defined on the patches.
template <class... Patches>
void initialise_all_functions(MultipatchField<DFieldOnPatch_host, Patches...> const& functions)
{
    (initialise_2D_function<typename Patches::Grid1, typename Patches::Grid2>(
             functions.template get<Patches>()),
     ...);
}


struct SingleInterfaceDerivativesCalculatorCollectionTest : public ::testing::Test
{
    // DEFINE BOUNDARIES OF THE DOMAINS ----------------------------------------------------------
    // patches 1 ---------------------------------
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


    // SingleInterfaceDerivativesCalculators for interfaces along y.
    SingleInterfaceDerivativesCalculator<Interface_12> const derivatives_calculator_1_2;
    SingleInterfaceDerivativesCalculator<Interface_23> const derivatives_calculator_2_3;

public:
    SingleInterfaceDerivativesCalculatorCollectionTest()
        : idx_range_x1(SplineInterpPointsX<1, BCG, BCH>::template get_domain<GridX<1>>())
        , idx_range_x2(SplineInterpPointsX<2, BCH, BCH>::template get_domain<GridX<2>>())
        , idx_range_x3(SplineInterpPointsX<3, BCH, BCG>::template get_domain<GridX<3>>())
        , idx_range_y1(SplineInterpPointsY<1, BCG, BCH>::template get_domain<GridY<1>>())
        , idx_range_y2(SplineInterpPointsY<2, BCH, BCH>::template get_domain<GridY<2>>())
        , idx_range_y3(SplineInterpPointsY<3, BCH, BCG>::template get_domain<GridY<3>>())
        , idx_range_xy1(idx_range_x1, idx_range_y1)
        , idx_range_xy2(idx_range_x2, idx_range_y2)
        , idx_range_xy3(idx_range_x3, idx_range_y3)
        , derivatives_calculator_1_2(idx_range_xy1, idx_range_xy2, ddc::BoundCond::GREVILLE)
        , derivatives_calculator_2_3(
                  idx_range_xy2,
                  idx_range_xy3,
                  ddc::BoundCond::HERMITE,
                  ddc::BoundCond::GREVILLE)
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

        std::vector<Coord<X<1>>> interpolation_points_x1
                = get_interpolation_points_add_one_on_left(break_points_x1);
        std::vector<Coord<X<2>>> interpolation_points_x2 = break_points_x2;
        std::vector<Coord<X<3>>> interpolation_points_x3
                = get_interpolation_points_add_one_on_right(break_points_x3);

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
    }



    // TEST OPERATORS ============================================================================
    template <class DerivativesCalculatorType>
    void check_get_coeff_deriv_patch_1_2_call_single(
            DerivativesCalculatorType const& deriv_calculators_tested,
            DerivativesCalculatorType const& deriv_calculators_expected)
    {
        EXPECT_EQ(
                deriv_calculators_tested.get_coeff_deriv_patch_1(),
                deriv_calculators_expected.get_coeff_deriv_patch_1());
        EXPECT_EQ(
                deriv_calculators_tested.get_coeff_deriv_patch_2(),
                deriv_calculators_expected.get_coeff_deriv_patch_2());
    }


    template <class DerivativesCalculatorType, class... Patches>
    void check_get_function_coefficients_call_single(
            DerivativesCalculatorType const& deriv_calculators_tested,
            DerivativesCalculatorType const& deriv_calculators_expected,
            MultipatchField<DFieldOnPatch_host, Patches...> const& functions)
    {
        using InterfaceI = typename DerivativesCalculatorType::associated_interface;
        using PatchL = typename InterfaceI::Edge1::associated_patch;
        using PatchR = typename InterfaceI::Edge2::associated_patch;

        using GridParL = typename InterfaceI::Edge1::parallel_grid;
        using GridParR = typename InterfaceI::Edge2::parallel_grid;

        DFieldOnPatch_host<PatchL> function_L = functions.template get<PatchL>();
        DFieldOnPatch_host<PatchR> function_R = functions.template get<PatchR>();

        // Pick a random index between 0 and 10.
        Idx<GridParL> idx_slice_L(3);
        Idx<GridParR> idx_slice_R(3);

        EXPECT_EQ(
                deriv_calculators_tested.get_function_coefficients(
                        get_const_field(function_L[idx_slice_L]),
                        get_const_field(function_R[idx_slice_R])),
                deriv_calculators_expected.get_function_coefficients(
                        get_const_field(function_L[idx_slice_L]),
                        get_const_field(function_R[idx_slice_R])));
    }



    template <class DerivativesCalculatorCollection, class... DerivativesCalculatorType>
    void check_function_call(
            MultipatchField<DFieldOnPatch_host, Patch1, Patch2, Patch3> const& functions,
            DerivativesCalculatorCollection const& deriv_calculators_collection,
            DerivativesCalculatorType const&... deriv_calculators)
    {
        std::tuple<DerivativesCalculatorType const&...> deriv_calculators_tuple(
                deriv_calculators...);

        (check_get_coeff_deriv_patch_1_2_call_single(
                 deriv_calculators_collection
                         .template get<typename DerivativesCalculatorType::associated_interface>(),
                 std::get<DerivativesCalculatorType const&>(deriv_calculators_tuple)),
         ...);


        (check_get_function_coefficients_call_single(
                 deriv_calculators_collection
                         .template get<typename DerivativesCalculatorType::associated_interface>(),
                 std::get<DerivativesCalculatorType const&>(deriv_calculators_tuple),
                 functions),
         ...);
    }
};

} // end namespace


TEST_F(SingleInterfaceDerivativesCalculatorCollectionTest, CheckCallToOperators)
{
    // Order in sequences ------------------------------------------------------------------------
    SingleInterfaceDerivativesCalculatorCollection
            deriv_calculators_collect(derivatives_calculator_1_2, derivatives_calculator_2_3);

    // Instantiate test function values ==========================================================
    host_t<DFieldMemOnPatch<Patch1>> function_1_alloc(idx_range_xy1);
    host_t<DFieldMemOnPatch<Patch2>> function_2_alloc(idx_range_xy2);
    host_t<DFieldMemOnPatch<Patch3>> function_3_alloc(idx_range_xy3);

    host_t<DFieldOnPatch<Patch1>> function_1(get_field(function_1_alloc));
    host_t<DFieldOnPatch<Patch2>> function_2(get_field(function_2_alloc));
    host_t<DFieldOnPatch<Patch3>> function_3(get_field(function_3_alloc));

    // Collect the function values.
    MultipatchField<DFieldOnPatch_host, Patch1, Patch2, Patch3>
            functions(function_1, function_2, function_3);

    // Initialise the function values
    initialise_all_functions(functions);

    // Test SingleInterfaceDerivativesCalculatorCollection =======================================
    check_function_call(
            functions,
            deriv_calculators_collect,
            derivatives_calculator_1_2,
            derivatives_calculator_2_3);
}