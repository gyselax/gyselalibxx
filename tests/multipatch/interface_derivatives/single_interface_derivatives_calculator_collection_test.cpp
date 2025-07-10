// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "9patches_2d_periodic_strips_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "interface.hpp"
#include "interface_derivative_matrix.hpp"
#include "interface_derivatives_test_utils.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "single_interface_derivatives_calculator.hpp"
#include "single_interface_derivatives_calculator_collection.hpp"
#include "view.hpp"


/*
    Test SingleInterfaceDerivativesCalculatorCollection on the following geometry:

        |  1  |  2  |  3  |  1 ...
        -------------------
        |  4  |  5  |  6  |  4 ...
        -------------------
        |  7  |  8  |  9  |  7 ... 

        with the global X dimension periodic and the global Y spline
        with additional points as closure condition (ddc::BoundCond::GREVILLE).
*/

namespace {
// Multi-patch tags ---
using namespace periodic_strips_non_uniform_2d_9patches;

using HostExecSpace = Kokkos::DefaultHostExecutionSpace;


// Interpolation points type for the patches.
template <std::size_t PatchIdx>
using SplineInterpPointsX = ddcHelper::NonUniformInterpolationPoints<
        BSplinesX<PatchIdx>,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

template <std::size_t PatchIdx, ddc::BoundCond BoundCondMin, ddc::BoundCond BoundCondMax>
using SplineInterpPointsY
        = ddcHelper::NonUniformInterpolationPoints<BSplinesY<1>, BoundCondMin, BoundCondMax>;

// USEFUL FUNCTIONS ==============================================================================
/// @brief Initialise the function with f(x,y) = cos(2/3*pi*x)sin(y -0.3).
template <class Grid1, class Grid2>
void initialise_2D_function(host_t<DField<IdxRange<Grid1, Grid2>>> function)
{
    ddc::for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
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
    // patches 1 | 4 | 7  dim X ------------------
    static constexpr Coord<X> x1_min = Coord<X>(0.0);
    static constexpr Coord<X> x1_max = Coord<X>(1.0);
    static constexpr IdxStep<GridX<1>> x1_ncells = IdxStep<GridX<1>>(10);

    static constexpr Coord<X> x4_min = x1_min;
    static constexpr Coord<X> x4_max = x1_max;
    static constexpr IdxStep<GridX<4>> x4_ncells = IdxStep<GridX<4>>(10);

    static constexpr Coord<X> x7_min = x1_min;
    static constexpr Coord<X> x7_max = x1_max;
    static constexpr IdxStep<GridX<7>> x7_ncells = IdxStep<GridX<7>>(10);

    // patches 2 | 5 | 8  dim X ------------------
    static constexpr Coord<X> x2_min = Coord<X>(1.0);
    static constexpr Coord<X> x2_max = Coord<X>(2.0);
    static constexpr IdxStep<GridX<2>> x2_ncells = IdxStep<GridX<2>>(10);

    static constexpr Coord<X> x5_min = x2_min;
    static constexpr Coord<X> x5_max = x2_max;
    static constexpr IdxStep<GridX<5>> x5_ncells = IdxStep<GridX<5>>(10);

    static constexpr Coord<X> x8_min = x2_min;
    static constexpr Coord<X> x8_max = x2_max;
    static constexpr IdxStep<GridX<8>> x8_ncells = IdxStep<GridX<8>>(10);

    // patches 3 | 6 | 9  dim X ------------------
    static constexpr Coord<X> x3_min = Coord<X>(2.0);
    static constexpr Coord<X> x3_max = Coord<X>(3.0);
    static constexpr IdxStep<GridX<3>> x3_ncells = IdxStep<GridX<3>>(10);

    static constexpr Coord<X> x6_min = x3_min;
    static constexpr Coord<X> x6_max = x3_max;
    static constexpr IdxStep<GridX<6>> x6_ncells = IdxStep<GridX<6>>(10);

    static constexpr Coord<X> x9_min = x3_min;
    static constexpr Coord<X> x9_max = x3_max;
    static constexpr IdxStep<GridX<9>> x9_ncells = IdxStep<GridX<9>>(10);

    // patches 1 | 2 | 3  dim Y --------------------
    static constexpr Coord<Y> y1_min = Coord<Y>(2.0);
    static constexpr Coord<Y> y1_max = Coord<Y>(3.0);
    static constexpr IdxStep<GridY<1>> y1_ncells = IdxStep<GridY<1>>(10);

    static constexpr Coord<Y> y2_min = y1_min;
    static constexpr Coord<Y> y2_max = y1_max;
    static constexpr IdxStep<GridY<2>> y2_ncells = IdxStep<GridY<2>>(10);

    static constexpr Coord<Y> y3_min = y1_min;
    static constexpr Coord<Y> y3_max = y1_max;
    static constexpr IdxStep<GridY<3>> y3_ncells = IdxStep<GridY<3>>(10);

    // patches 4 | 5 | 6  dim Y --------------------
    static constexpr Coord<Y> y4_min = Coord<Y>(1.0);
    static constexpr Coord<Y> y4_max = Coord<Y>(2.0);
    static constexpr IdxStep<GridY<4>> y4_ncells = IdxStep<GridY<4>>(10);

    static constexpr Coord<Y> y5_min = y4_min;
    static constexpr Coord<Y> y5_max = y4_max;
    static constexpr IdxStep<GridY<5>> y5_ncells = IdxStep<GridY<5>>(10);

    static constexpr Coord<Y> y6_min = y4_min;
    static constexpr Coord<Y> y6_max = y4_max;
    static constexpr IdxStep<GridY<6>> y6_ncells = IdxStep<GridY<6>>(10);

    // patches 7 | 8 | 9  dim Y --------------------
    static constexpr Coord<Y> y7_min = Coord<Y>(0.0);
    static constexpr Coord<Y> y7_max = Coord<Y>(1.0);
    static constexpr IdxStep<GridY<7>> y7_ncells = IdxStep<GridY<7>>(10);

    static constexpr Coord<Y> y8_min = y7_min;
    static constexpr Coord<Y> y8_max = y7_max;
    static constexpr IdxStep<GridY<8>> y8_ncells = IdxStep<GridY<8>>(10);

    static constexpr Coord<Y> y9_min = y7_min;
    static constexpr Coord<Y> y9_max = y7_max;
    static constexpr IdxStep<GridY<9>> y9_ncells = IdxStep<GridY<9>>(10);

    static constexpr ddc::BoundCond BcH = ddc::BoundCond::HERMITE;
    static constexpr ddc::BoundCond BcG = ddc::BoundCond::GREVILLE;

protected:
    const IdxRange<GridX<1>> idx_range_x1;
    const IdxRange<GridX<2>> idx_range_x2;
    const IdxRange<GridX<3>> idx_range_x3;
    const IdxRange<GridX<4>> idx_range_x4;
    const IdxRange<GridX<5>> idx_range_x5;
    const IdxRange<GridX<6>> idx_range_x6;
    const IdxRange<GridX<7>> idx_range_x7;
    const IdxRange<GridX<8>> idx_range_x8;
    const IdxRange<GridX<9>> idx_range_x9;


    const IdxRange<GridY<1>> idx_range_y1;
    const IdxRange<GridY<2>> idx_range_y2;
    const IdxRange<GridY<3>> idx_range_y3;
    const IdxRange<GridY<4>> idx_range_y4;
    const IdxRange<GridY<5>> idx_range_y5;
    const IdxRange<GridY<6>> idx_range_y6;
    const IdxRange<GridY<7>> idx_range_y7;
    const IdxRange<GridY<8>> idx_range_y8;
    const IdxRange<GridY<9>> idx_range_y9;


    const IdxRange<GridX<1>, GridY<1>> idx_range_xy1;
    const IdxRange<GridX<2>, GridY<2>> idx_range_xy2;
    const IdxRange<GridX<3>, GridY<3>> idx_range_xy3;
    const IdxRange<GridX<4>, GridY<4>> idx_range_xy4;
    const IdxRange<GridX<5>, GridY<5>> idx_range_xy5;
    const IdxRange<GridX<6>, GridY<6>> idx_range_xy6;
    const IdxRange<GridX<7>, GridY<7>> idx_range_xy7;
    const IdxRange<GridX<8>, GridY<8>> idx_range_xy8;
    const IdxRange<GridX<9>, GridY<9>> idx_range_xy9;


    // SingleInterfaceDerivativesCalculators for interfaces along y (periodic).
    SingleInterfaceDerivativesCalculator<Interface_1_2> const derivatives_calculator_1_2;
    SingleInterfaceDerivativesCalculator<Interface_2_3> const derivatives_calculator_2_3;
    SingleInterfaceDerivativesCalculator<Interface_3_1> const derivatives_calculator_3_1;

    SingleInterfaceDerivativesCalculator<Interface_4_5> const derivatives_calculator_4_5;
    SingleInterfaceDerivativesCalculator<Interface_5_6> const derivatives_calculator_5_6;
    SingleInterfaceDerivativesCalculator<Interface_6_4> const derivatives_calculator_6_4;

    SingleInterfaceDerivativesCalculator<Interface_7_8> const derivatives_calculator_7_8;
    SingleInterfaceDerivativesCalculator<Interface_8_9> const derivatives_calculator_8_9;
    SingleInterfaceDerivativesCalculator<Interface_9_7> const derivatives_calculator_9_7;

    // SingleInterfaceDerivativesCalculators for interfaces along x.
    SingleInterfaceDerivativesCalculator<
            Interface_1_4,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::HERMITE> const derivatives_calculator_1_4;
    SingleInterfaceDerivativesCalculator<
            Interface_4_7,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::GREVILLE> const derivatives_calculator_4_7;

    SingleInterfaceDerivativesCalculator<
            Interface_2_5,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::HERMITE> const derivatives_calculator_2_5;
    SingleInterfaceDerivativesCalculator<
            Interface_5_8,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::GREVILLE> const derivatives_calculator_5_8;

    SingleInterfaceDerivativesCalculator<
            Interface_3_6,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::HERMITE> const derivatives_calculator_3_6;
    SingleInterfaceDerivativesCalculator<
            Interface_6_9,
            ddc::BoundCond::HERMITE,
            ddc::BoundCond::GREVILLE> const derivatives_calculator_6_9;

public:
    SingleInterfaceDerivativesCalculatorCollectionTest()
        : idx_range_x1(SplineInterpPointsX<1>::template get_domain<GridX<1>>())
        , idx_range_x2(SplineInterpPointsX<2>::template get_domain<GridX<2>>())
        , idx_range_x3(SplineInterpPointsX<3>::template get_domain<GridX<3>>())
        , idx_range_x4(SplineInterpPointsX<4>::template get_domain<GridX<4>>())
        , idx_range_x5(SplineInterpPointsX<5>::template get_domain<GridX<5>>())
        , idx_range_x6(SplineInterpPointsX<6>::template get_domain<GridX<6>>())
        , idx_range_x7(SplineInterpPointsX<7>::template get_domain<GridX<7>>())
        , idx_range_x8(SplineInterpPointsX<8>::template get_domain<GridX<8>>())
        , idx_range_x9(SplineInterpPointsX<9>::template get_domain<GridX<9>>())
        , idx_range_y1(SplineInterpPointsY<1, BcH, BcG>::template get_domain<GridY<1>>())
        , idx_range_y2(SplineInterpPointsY<2, BcH, BcG>::template get_domain<GridY<2>>())
        , idx_range_y3(SplineInterpPointsY<3, BcH, BcG>::template get_domain<GridY<3>>())
        , idx_range_y4(SplineInterpPointsY<4, BcH, BcH>::template get_domain<GridY<4>>())
        , idx_range_y5(SplineInterpPointsY<5, BcH, BcH>::template get_domain<GridY<5>>())
        , idx_range_y6(SplineInterpPointsY<6, BcH, BcH>::template get_domain<GridY<6>>())
        , idx_range_y7(SplineInterpPointsY<7, BcG, BcH>::template get_domain<GridY<7>>())
        , idx_range_y8(SplineInterpPointsY<8, BcG, BcH>::template get_domain<GridY<8>>())
        , idx_range_y9(SplineInterpPointsY<9, BcG, BcH>::template get_domain<GridY<9>>())
        , idx_range_xy1(idx_range_x1, idx_range_y1)
        , idx_range_xy2(idx_range_x2, idx_range_y2)
        , idx_range_xy3(idx_range_x3, idx_range_y3)
        , idx_range_xy4(idx_range_x4, idx_range_y4)
        , idx_range_xy5(idx_range_x5, idx_range_y5)
        , idx_range_xy6(idx_range_x6, idx_range_y6)
        , idx_range_xy7(idx_range_x7, idx_range_y7)
        , idx_range_xy8(idx_range_x8, idx_range_y8)
        , idx_range_xy9(idx_range_x9, idx_range_y9)
        , derivatives_calculator_1_2(idx_range_xy1, idx_range_xy2)
        , derivatives_calculator_2_3(idx_range_xy2, idx_range_xy3)
        , derivatives_calculator_3_1(idx_range_xy3, idx_range_xy1)
        , derivatives_calculator_4_5(idx_range_xy4, idx_range_xy5)
        , derivatives_calculator_5_6(idx_range_xy5, idx_range_xy6)
        , derivatives_calculator_6_4(idx_range_xy6, idx_range_xy4)
        , derivatives_calculator_7_8(idx_range_xy7, idx_range_xy8)
        , derivatives_calculator_8_9(idx_range_xy8, idx_range_xy9)
        , derivatives_calculator_9_7(idx_range_xy9, idx_range_xy7)
        , derivatives_calculator_1_4(idx_range_xy1, idx_range_xy4)
        , derivatives_calculator_4_7(idx_range_xy4, idx_range_xy7)
        , derivatives_calculator_2_5(idx_range_xy2, idx_range_xy5)
        , derivatives_calculator_5_8(idx_range_xy5, idx_range_xy8)
        , derivatives_calculator_3_6(idx_range_xy3, idx_range_xy6)
        , derivatives_calculator_6_9(idx_range_xy6, idx_range_xy9)
    {
    }


    // INITIALISE DOMAINS ------------------------------------------------------------------------
    static void SetUpTestSuite()
    {
        // Creating of meshes and supports .......................................................
        std::vector<Coord<X>> break_points_x147
                = build_random_non_uniform_break_points(x1_min, x1_max, x1_ncells);
        std::vector<Coord<X>> break_points_x258
                = build_random_non_uniform_break_points(x2_min, x2_max, x2_ncells);
        std::vector<Coord<X>> break_points_x369
                = build_random_non_uniform_break_points(x3_min, x3_max, x3_ncells);

        std::vector<Coord<Y>> break_points_y123
                = build_random_non_uniform_break_points(y1_min, y1_max, y1_ncells);
        std::vector<Coord<Y>> break_points_y456
                = build_random_non_uniform_break_points(y4_min, y4_max, y4_ncells);
        std::vector<Coord<Y>> break_points_y789
                = build_random_non_uniform_break_points(y7_min, y7_max, y7_ncells);

        std::vector<Coord<X>> interpolation_points_x147 = break_points_x147;
        std::vector<Coord<X>> interpolation_points_x258 = break_points_x258;
        std::vector<Coord<X>> interpolation_points_x369 = break_points_x369;

        std::vector<Coord<Y>> interpolation_points_y123
                = get_interpolation_points_add_one_on_right(break_points_y123);
        std::vector<Coord<Y>> interpolation_points_y456 = break_points_y456;
        std::vector<Coord<Y>> interpolation_points_y789
                = get_interpolation_points_add_one_on_left(break_points_y789);


        // Patch 1 ...............................................................................
        ddc::init_discrete_space<BSplinesX<1>>(break_points_x147);
        ddc::init_discrete_space<BSplinesY<1>>(break_points_y123);

        ddc::init_discrete_space<GridX<1>>(interpolation_points_x147);
        ddc::init_discrete_space<GridY<1>>(interpolation_points_y123);

        // Patch 2 ...............................................................................
        ddc::init_discrete_space<BSplinesX<2>>(break_points_x258);
        ddc::init_discrete_space<BSplinesY<2>>(break_points_y123);

        ddc::init_discrete_space<GridX<2>>(interpolation_points_x258);
        ddc::init_discrete_space<GridY<2>>(interpolation_points_y123);

        // Patch 3 ...............................................................................
        ddc::init_discrete_space<BSplinesX<3>>(break_points_x369);
        ddc::init_discrete_space<BSplinesY<3>>(break_points_y123);

        ddc::init_discrete_space<GridX<3>>(interpolation_points_x369);
        ddc::init_discrete_space<GridY<3>>(interpolation_points_y123);

        // Patch 4 ...............................................................................
        ddc::init_discrete_space<BSplinesX<4>>(break_points_x147);
        ddc::init_discrete_space<BSplinesY<4>>(break_points_y456);

        ddc::init_discrete_space<GridX<4>>(interpolation_points_x147);
        ddc::init_discrete_space<GridY<4>>(interpolation_points_y456);

        // Patch 5 ...............................................................................
        ddc::init_discrete_space<BSplinesX<5>>(break_points_x258);
        ddc::init_discrete_space<BSplinesY<5>>(break_points_y456);

        ddc::init_discrete_space<GridX<5>>(interpolation_points_x258);
        ddc::init_discrete_space<GridY<5>>(interpolation_points_y456);

        // Patch 6 ...............................................................................
        ddc::init_discrete_space<BSplinesX<6>>(break_points_x369);
        ddc::init_discrete_space<BSplinesY<6>>(break_points_y456);

        ddc::init_discrete_space<GridX<6>>(interpolation_points_x369);
        ddc::init_discrete_space<GridY<6>>(interpolation_points_y456);

        // Patch 7 ...............................................................................
        ddc::init_discrete_space<BSplinesX<7>>(break_points_x147);
        ddc::init_discrete_space<BSplinesY<7>>(break_points_y789);

        ddc::init_discrete_space<GridX<7>>(interpolation_points_x147);
        ddc::init_discrete_space<GridY<7>>(interpolation_points_y789);

        // Patch 8 ...............................................................................
        ddc::init_discrete_space<BSplinesX<8>>(break_points_x258);
        ddc::init_discrete_space<BSplinesY<8>>(break_points_y789);

        ddc::init_discrete_space<GridX<8>>(interpolation_points_x258);
        ddc::init_discrete_space<GridY<8>>(interpolation_points_y789);

        // Patch 9 ...............................................................................
        ddc::init_discrete_space<BSplinesX<9>>(break_points_x369);
        ddc::init_discrete_space<BSplinesY<9>>(break_points_y789);

        ddc::init_discrete_space<GridX<9>>(interpolation_points_x369);
        ddc::init_discrete_space<GridY<9>>(interpolation_points_y789);
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


    template <class DerivativesCalculatorType>
    void check_get_function_coefficients_call_single(
            DerivativesCalculatorType const& deriv_calculators_tested,
            DerivativesCalculatorType const& deriv_calculators_expected,
            MultipatchField<
                    DFieldOnPatch_host,
                    Patch1,
                    Patch2,
                    Patch3,
                    Patch4,
                    Patch5,
                    Patch6,
                    Patch7,
                    Patch8,
                    Patch9> const& functions)
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
            MultipatchField<
                    DFieldOnPatch_host,
                    Patch1,
                    Patch2,
                    Patch3,
                    Patch4,
                    Patch5,
                    Patch6,
                    Patch7,
                    Patch8,
                    Patch9> const& functions,
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


// TEST_F(SingleInterfaceDerivativesCalculatorCollectionTest, CheckThrowWrongInterface)
// {
//     SingleInterfaceDerivativesCalculatorCollection deriv_calculators_collect_123(
//             derivatives_calculator_1_2,
//             derivatives_calculator_2_3,
//             derivatives_calculator_3_1);

//     EXPECT_THROW(auto deriv_calculator = deriv_calculators_collect_123.template get<Interface_1_2>(), std::invalid_argument);
// }

TEST_F(SingleInterfaceDerivativesCalculatorCollectionTest, CheckCallToOperators)
{
    // Order in sequences ------------------------------------------------------------------------
    SingleInterfaceDerivativesCalculatorCollection deriv_calculators_collect_123(
            derivatives_calculator_1_2,
            derivatives_calculator_2_3,
            derivatives_calculator_3_1);

    SingleInterfaceDerivativesCalculatorCollection deriv_calculators_collect_456(
            derivatives_calculator_4_5,
            derivatives_calculator_5_6,
            derivatives_calculator_6_4);

    SingleInterfaceDerivativesCalculatorCollection deriv_calculators_collect_789(
            derivatives_calculator_7_8,
            derivatives_calculator_8_9,
            derivatives_calculator_9_7);

    SingleInterfaceDerivativesCalculatorCollection
            deriv_calculators_collect_147(derivatives_calculator_1_4, derivatives_calculator_4_7);

    SingleInterfaceDerivativesCalculatorCollection
            deriv_calculators_collect_258(derivatives_calculator_2_5, derivatives_calculator_5_8);

    SingleInterfaceDerivativesCalculatorCollection
            deriv_calculators_collect_369(derivatives_calculator_3_6, derivatives_calculator_6_9);


    // Instantiate test function values ==========================================================
    host_t<DFieldMemOnPatch<Patch1>> function_1_alloc(idx_range_xy1);
    host_t<DFieldMemOnPatch<Patch2>> function_2_alloc(idx_range_xy2);
    host_t<DFieldMemOnPatch<Patch3>> function_3_alloc(idx_range_xy3);
    host_t<DFieldMemOnPatch<Patch4>> function_4_alloc(idx_range_xy4);
    host_t<DFieldMemOnPatch<Patch5>> function_5_alloc(idx_range_xy5);
    host_t<DFieldMemOnPatch<Patch6>> function_6_alloc(idx_range_xy6);
    host_t<DFieldMemOnPatch<Patch7>> function_7_alloc(idx_range_xy7);
    host_t<DFieldMemOnPatch<Patch8>> function_8_alloc(idx_range_xy8);
    host_t<DFieldMemOnPatch<Patch9>> function_9_alloc(idx_range_xy9);

    host_t<DFieldOnPatch<Patch1>> function_1(get_field(function_1_alloc));
    host_t<DFieldOnPatch<Patch2>> function_2(get_field(function_2_alloc));
    host_t<DFieldOnPatch<Patch3>> function_3(get_field(function_3_alloc));
    host_t<DFieldOnPatch<Patch4>> function_4(get_field(function_4_alloc));
    host_t<DFieldOnPatch<Patch5>> function_5(get_field(function_5_alloc));
    host_t<DFieldOnPatch<Patch6>> function_6(get_field(function_6_alloc));
    host_t<DFieldOnPatch<Patch7>> function_7(get_field(function_7_alloc));
    host_t<DFieldOnPatch<Patch8>> function_8(get_field(function_8_alloc));
    host_t<DFieldOnPatch<Patch9>> function_9(get_field(function_9_alloc));

    // Collect the function values.
    MultipatchField<
            DFieldOnPatch_host,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            functions(
                    function_1,
                    function_2,
                    function_3,
                    function_4,
                    function_5,
                    function_6,
                    function_7,
                    function_8,
                    function_9);

    // Initialise the function values
    initialise_all_functions(functions);

    // Test SingleInterfaceDerivativesCalculatorCollection =======================================
    check_function_call(
            functions,
            deriv_calculators_collect_123,
            derivatives_calculator_1_2,
            derivatives_calculator_2_3,
            derivatives_calculator_3_1);

    check_function_call(
            functions,
            deriv_calculators_collect_456,
            derivatives_calculator_4_5,
            derivatives_calculator_5_6,
            derivatives_calculator_6_4);


    check_function_call(
            functions,
            deriv_calculators_collect_789,
            derivatives_calculator_7_8,
            derivatives_calculator_8_9,
            derivatives_calculator_9_7);

    check_function_call(
            functions,
            deriv_calculators_collect_147,
            derivatives_calculator_1_4,
            derivatives_calculator_4_7);

    check_function_call(
            functions,
            deriv_calculators_collect_258,
            derivatives_calculator_2_5,
            derivatives_calculator_5_8);

    check_function_call(
            functions,
            deriv_calculators_collect_369,
            derivatives_calculator_3_6,
            derivatives_calculator_6_9);
}
