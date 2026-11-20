// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "../periodic_strips_non_uniform_9patches_testing_tools.hpp"

#include "9patches_2d_periodic_strips_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "interface.hpp"
#include "interface_derivative_coefficients.hpp"
#include "interface_derivatives_test_utils.hpp"
#include "interfaces_derivative_calculator.hpp"
#include "interfaces_derivative_calculator_test_utils.hpp"
#include "interfaces_exact_derivative_calculator.hpp"
#include "linear_coord_transform.hpp"
#include "mesh_builder.hpp"
#include "multipatch_field_mem.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "orthogonal_coord_transforms.hpp"
#include "test_utils.hpp"
#include "types.hpp"
#include "view.hpp"


/*
    Test InterfacesDerivativeCalculator or InterfacesExactDerivativeCalculator on the following geometry:

        |  1  |  2  |  3  |  1 ...
        -------------------
        |  4  |  5  |  6  |  4 ...
        -------------------
        |  7  |  8  |  9  |  7 ... 

        with the global X dimension periodic and the global Y spline
        with additional points as closure condition (ddc::SplineBuilderClosure::GREVILLE).

    > test ddc::SplineBuilderClosure::PERIODIC boundary conditions. 
    > test ddc::SplineBuilderClosure::GREVILLE boundary conditions. 
    > test application on the X and Y directions. 
    > test application to compute first derivatives and cross-derivatives. 
    > test on a non-uniform patches. 
    > test with approximation formulation in InterfacesDerivCoeffs. 
    > test agreement between computed and global spline derivatives. 
    > test agreement between local and global splines.
*/

namespace {

template <class is_approximation_literal>
struct InterfacesDerivativeCalculatorGrevillePeriodicTest
    : public Multipatch9PatchesPeriodicStripsNonUniformTest<is_approximation_literal>
{
};

using is_approximation = std::tuple<std::true_type, std::false_type>;
using TestTypes = tuple_to_types_t<is_approximation>;
TYPED_TEST_SUITE(InterfacesDerivativeCalculatorGrevillePeriodicTest, TestTypes);

} // end namespace



TYPED_TEST(InterfacesDerivativeCalculatorGrevillePeriodicTest, CheckForPeriodicAndGrevilleBC)
{
    std::tuple coord_transforms {
            CoordTransformGroup<1>(),
            CoordTransformGroup<2>(),
            CoordTransformGroup<3>(),
            CoordTransformGroup<4>(),
            CoordTransformGroup<5>(),
            CoordTransformGroup<6>(),
            CoordTransformGroup<7>(),
            CoordTransformGroup<8>(),
            CoordTransformGroup<9>()};

    // Instantiate the derivatives calculators ---------------------------------------------------
    // InterfaceDerivCoeffs for interfaces along y (periodic).
    InterfaceDerivCoeffs<Interface_1_2> const
            deriv_coeffs_1_2(TestFixture::idx_range_xy1, TestFixture::idx_range_xy2);
    InterfaceDerivCoeffs<Interface_2_3> const
            deriv_coeffs_2_3(TestFixture::idx_range_xy2, TestFixture::idx_range_xy3);
    InterfaceDerivCoeffs<Interface_3_1> const
            deriv_coeffs_3_1(TestFixture::idx_range_xy3, TestFixture::idx_range_xy1);

    InterfaceDerivCoeffs<Interface_4_5> const
            deriv_coeffs_4_5(TestFixture::idx_range_xy4, TestFixture::idx_range_xy5);
    InterfaceDerivCoeffs<Interface_5_6> const
            deriv_coeffs_5_6(TestFixture::idx_range_xy5, TestFixture::idx_range_xy6);
    InterfaceDerivCoeffs<Interface_6_4> const
            deriv_coeffs_6_4(TestFixture::idx_range_xy6, TestFixture::idx_range_xy4);

    InterfaceDerivCoeffs<Interface_7_8> const
            deriv_coeffs_7_8(TestFixture::idx_range_xy7, TestFixture::idx_range_xy8);
    InterfaceDerivCoeffs<Interface_8_9> const
            deriv_coeffs_8_9(TestFixture::idx_range_xy8, TestFixture::idx_range_xy9);
    InterfaceDerivCoeffs<Interface_9_7> const
            deriv_coeffs_9_7(TestFixture::idx_range_xy9, TestFixture::idx_range_xy7);

    // InterfaceDerivCoeffs for interfaces along x.
    constexpr std::size_t nb_chosen_cells = 9;

    constexpr bool is_approximation = TestFixture::is_approximation;

    // InterfaceDerivCoeffs for interfaces along x.
    InterfaceDerivCoeffs<Interface_1_4> const deriv_coeffs_1_4
            = build_interface<is_approximation, Interface_1_4>(
                    TestFixture::idx_range_xy1,
                    TestFixture::idx_range_xy4,
                    nb_chosen_cells,
                    ddc::SplineBuilderClosure::GREVILLE,
                    ddc::SplineBuilderClosure::HERMITE);
    InterfaceDerivCoeffs<Interface_4_7> const deriv_coeffs_4_7
            = build_interface<is_approximation, Interface_4_7>(
                    TestFixture::idx_range_xy4,
                    TestFixture::idx_range_xy7,
                    nb_chosen_cells,
                    ddc::SplineBuilderClosure::HERMITE,
                    ddc::SplineBuilderClosure::GREVILLE);

    InterfaceDerivCoeffs<Interface_2_5> const deriv_coeffs_2_5
            = build_interface<is_approximation, Interface_2_5>(
                    TestFixture::idx_range_xy2,
                    TestFixture::idx_range_xy5,
                    nb_chosen_cells,
                    ddc::SplineBuilderClosure::GREVILLE,
                    ddc::SplineBuilderClosure::HERMITE);
    InterfaceDerivCoeffs<Interface_5_8> const deriv_coeffs_5_8
            = build_interface<is_approximation, Interface_5_8>(
                    TestFixture::idx_range_xy5,
                    TestFixture::idx_range_xy8,
                    nb_chosen_cells,
                    ddc::SplineBuilderClosure::HERMITE,
                    ddc::SplineBuilderClosure::GREVILLE);

    InterfaceDerivCoeffs<Interface_3_6> const deriv_coeffs_3_6
            = build_interface<is_approximation, Interface_3_6>(
                    TestFixture::idx_range_xy3,
                    TestFixture::idx_range_xy6,
                    nb_chosen_cells,
                    ddc::SplineBuilderClosure::GREVILLE,
                    ddc::SplineBuilderClosure::HERMITE);
    InterfaceDerivCoeffs<Interface_6_9> const deriv_coeffs_6_9
            = build_interface<is_approximation, Interface_6_9>(
                    TestFixture::idx_range_xy6,
                    TestFixture::idx_range_xy9,
                    nb_chosen_cells,
                    ddc::SplineBuilderClosure::HERMITE,
                    ddc::SplineBuilderClosure::GREVILLE);

    // Collect the index ranges ------------------------------------------------------------------
    MultipatchType<
            IdxRangeOnPatch,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            idx_ranges(
                    TestFixture::idx_range_xy1,
                    TestFixture::idx_range_xy2,
                    TestFixture::idx_range_xy3,
                    TestFixture::idx_range_xy4,
                    TestFixture::idx_range_xy5,
                    TestFixture::idx_range_xy6,
                    TestFixture::idx_range_xy7,
                    TestFixture::idx_range_xy8,
                    TestFixture::idx_range_xy9);

    // Instantiate the derivatives calculators ---------------------------------------------------
    InterfacesDerivativeCalculatorType<
            is_approximation,
            Connectivity,
            GridX<1>,
            ddc::SplineBuilderClosure::PERIODIC>
            deriv_calculator_123(idx_ranges, deriv_coeffs_1_2, deriv_coeffs_2_3, deriv_coeffs_3_1);

    InterfacesDerivativeCalculatorType<
            is_approximation,
            Connectivity,
            GridX<4>,
            ddc::SplineBuilderClosure::PERIODIC>
            deriv_calculator_456(idx_ranges, deriv_coeffs_4_5, deriv_coeffs_5_6, deriv_coeffs_6_4);

    InterfacesDerivativeCalculatorType<
            is_approximation,
            Connectivity,
            GridX<7>,
            ddc::SplineBuilderClosure::PERIODIC>
            deriv_calculator_789(idx_ranges, deriv_coeffs_7_8, deriv_coeffs_8_9, deriv_coeffs_9_7);

    InterfacesDerivativeCalculatorType<
            is_approximation,
            Connectivity,
            GridY<1>,
            ddc::SplineBuilderClosure::GREVILLE>
            deriv_calculator_147(idx_ranges, deriv_coeffs_1_4, deriv_coeffs_4_7);

    InterfacesDerivativeCalculatorType<
            is_approximation,
            Connectivity,
            GridY<2>,
            ddc::SplineBuilderClosure::GREVILLE>
            deriv_calculator_258(idx_ranges, deriv_coeffs_2_5, deriv_coeffs_5_8);

    InterfacesDerivativeCalculatorType<
            is_approximation,
            Connectivity,
            GridY<3>,
            ddc::SplineBuilderClosure::GREVILLE>
            deriv_calculator_369(idx_ranges, deriv_coeffs_3_6, deriv_coeffs_6_9);

    // Instantiate DerivField ====================================================================
    MultipatchFieldMem<
            DerivFieldMemOnPatch_host,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            functions_and_derivs_alloc(idx_ranges);
    MultipatchField<
            DerivFieldOnPatch_host,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            functions_and_derivs(functions_and_derivs_alloc);

    // Instantiate the field of global function values. No derivatives needed.
    host_t<DFieldMem<IdxRange<GridXg, GridYg>>> function_g_alloc(TestFixture::idx_range_xy_g);
    host_t<DField<IdxRange<GridXg, GridYg>>> function_g = get_field(function_g_alloc);

    // Initialise the data =======================================================================
    // --- the function values.
    initialise_all_functions<Xg, Yg>(functions_and_derivs, coord_transforms);
    initialise_2D_function(function_g);

    // --- the first derivatives computed from the function values.
    deriv_calculator_123(Idx<ddc::Deriv<X<1>>>(1), functions_and_derivs);
    deriv_calculator_456(Idx<ddc::Deriv<X<4>>>(1), functions_and_derivs);
    deriv_calculator_789(Idx<ddc::Deriv<X<7>>>(1), functions_and_derivs);

    deriv_calculator_147(Idx<ddc::Deriv<Y<1>>>(1), functions_and_derivs);
    deriv_calculator_258(Idx<ddc::Deriv<Y<2>>>(1), functions_and_derivs);
    deriv_calculator_369(Idx<ddc::Deriv<Y<3>>>(1), functions_and_derivs);

    // --- the cross-derivatives computed from the first derivatives.
    deriv_calculator_147(Idx<ddc::Deriv<X<1>>, ddc::Deriv<Y<1>>>(1, 1), functions_and_derivs);
    deriv_calculator_258(Idx<ddc::Deriv<X<2>>, ddc::Deriv<Y<2>>>(1, 1), functions_and_derivs);
    deriv_calculator_369(Idx<ddc::Deriv<X<3>>, ddc::Deriv<Y<3>>>(1, 1), functions_and_derivs);

    // Test the values of the derivatives ========================================================
    // --- Define an equivalent global spline.
    // Build global spline representation ---
    GlobalSplineBuilder builder_g(TestFixture::idx_range_xy_g);

    host_t<DFieldMem<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef_alloc(
            builder_g.batched_spline_domain(TestFixture::idx_range_xy_g));
    host_t<DField<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef
            = get_field(function_g_coef_alloc);

    builder_g(function_g_coef, get_const_field(function_g));

    // Global spline evaluator ---
    ddc::ConstantExtrapolationRule<Yg> bc_ymin_g(TestFixture::yg_min);
    ddc::ConstantExtrapolationRule<Yg> bc_ymax_g(TestFixture::yg_max);
    ddc::PeriodicExtrapolationRule<Xg> bc_x_g;
    GlobalSplineEvaluator evaluator_g(bc_x_g, bc_x_g, bc_ymin_g, bc_ymax_g);

    // Check each derivatives ---
    using PatchSeqLowerBound = ddc::TypeSeq<Patch7, Patch8, Patch9>;
    using PatchSeqUpperBound = ddc::TypeSeq<Patch1, Patch2, Patch3>;

    const double TOL_deriv = is_approximation ? 1e-4 : 1e-11;
    const double TOL_cross_deriv = is_approximation ? 1e-4 : 1e-10;
    const double TOL_spline = is_approximation ? 1e-7 : 1e-12;

    check_all_x_derivatives(
            functions_and_derivs,
            evaluator_g,
            get_const_field(function_g_coef),
            coord_transforms,
            TOL_deriv);
    check_all_y_derivatives<PatchSeqLowerBound, PatchSeqUpperBound>(
            functions_and_derivs,
            evaluator_g,
            get_const_field(function_g_coef),
            coord_transforms,
            TOL_deriv);
    check_all_xy_derivatives<PatchSeqLowerBound, PatchSeqUpperBound>(
            functions_and_derivs,
            evaluator_g,
            get_const_field(function_g_coef),
            coord_transforms,
            TOL_cross_deriv);

    // Check the whole spline representations ---
    check_all_spline_representation_agreement<PatchSeqLowerBound, PatchSeqUpperBound>(
            functions_and_derivs,
            evaluator_g,
            get_const_field(function_g_coef),
            coord_transforms,
            TOL_spline);
}
