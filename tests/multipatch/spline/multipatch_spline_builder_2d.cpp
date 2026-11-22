// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "../interface_derivatives/interface_derivatives_test_utils.hpp"
#include "../interface_derivatives/interfaces_derivative_calculator_test_utils.hpp"
#include "../periodic_strips_non_uniform_9patches_testing_tools.hpp"

#include "9patches_2d_periodic_strips_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "interface.hpp"
#include "interface_derivative_coefficients.hpp"
#include "interfaces_derivative_calculator.hpp"
#include "linear_coord_transform.hpp"
#include "mesh_builder.hpp"
#include "multipatch_field_mem.hpp"
#include "multipatch_spline_builder_2d.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "orthogonal_coord_transforms.hpp"
#include "types.hpp"
#include "view.hpp"


/*
    Test MultipatchSplineBuilder2D on the following geometry:

        |  1  |  2  |  3  |  1 ...
        -------------------
        |  4  |  5  |  6  |  4 ...
        -------------------
        |  7  |  8  |  9  |  7 ... 

        with the global X dimension periodic and the global Y spline
        with additional points as closure condition (ddc::SplineBuilderClosure::GREVILLE).
*/
namespace {

struct MultipatchSplineBuilder2DTest
    : public Multipatch9PatchesPeriodicStripsNonUniformTest<std::false_type>
{
};

template <
        class P,
        ddc::SplineBuilderClosure BC1min,
        ddc::SplineBuilderClosure BC1max,
        ddc::SplineBuilderClosure BC2min,
        ddc::SplineBuilderClosure BC2max>
using LocalSplineBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplines1OnPatch<P>,
        BSplines2OnPatch<P>,
        Grid1OnPatch<P>,
        Grid2OnPatch<P>,
        BC1min,
        BC1max,
        BC2min,
        BC2max,
        ddc::SplineSolver::LAPACK>;

template <class P>
using LocalSplineEvaluator = ddc::SplineEvaluatorND<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        ddc::TypeSeq<BSplines1OnPatch<P>, BSplines2OnPatch<P>>,
        ddc::TypeSeq<Grid1OnPatch<P>, Grid2OnPatch<P>>,
        ddc::TypeSeq<
                ddc::NullExtrapolationRule,
                ddc::NullExtrapolationRule,
                ddc::NullExtrapolationRule,
                ddc::NullExtrapolationRule>>;


template <class P, class CoordTransformType>
void check_spline(
        std::vector<double>& random_values_x,
        std::vector<double>& random_values_y,
        SplineCoeffOnPatch_2D<P> const& spline_l,
        host_t<DField<IdxRange<BSplinesXg, BSplinesYg>>> const& spline_g,
        LocalSplineEvaluator<P> const& evaluator_l,
        GlobalSplineEvaluator const& evaluator_g,
        CoordTransformType const& coord_transform,
        double const TOL = 1e-14)
{
    inverse_mapping_t<typename CoordTransformType::XYTransform> coord_transform_l_to_g_2d
            = coord_transform.coord_transform.get_inverse_mapping();

    // Evaluate and compare the local and global spline representations ----------------------
    double xmin = ddc::discrete_space<typename P::BSplines1>().rmin();
    double ymin = ddc::discrete_space<typename P::BSplines2>().rmin();
    double length_x = ddc::discrete_space<typename P::BSplines1>().length();
    double length_y = ddc::discrete_space<typename P::BSplines2>().length();

    for (double const val_x : random_values_x) {
        for (double const val_y : random_values_y) {
            typename P::Coord12 eval_point_l(xmin + val_x * length_x, ymin + val_y * length_y);
            Coord<Xg, Yg> eval_point_g(coord_transform_l_to_g_2d(eval_point_l));
            double local_spline = evaluator_l(eval_point_l, get_const_field(spline_l));
            double global_spline = evaluator_g(eval_point_g, get_const_field(spline_g));
            EXPECT_NEAR(local_spline, global_spline, TOL);
        }
    }
};

template <class... Patches, class CoordTransformTuple, class SplineEvaluatorPatchTuple>
void check_all_splines(
        std::vector<double>& random_values_x,
        std::vector<double>& random_values_y,
        MultipatchField<SplineCoeffOnPatch_2D, Patches...> const& splines_l,
        host_t<DField<IdxRange<BSplinesXg, BSplinesYg>>> const& spline_g,
        SplineEvaluatorPatchTuple const& evaluators_l,
        GlobalSplineEvaluator const& evaluator_g,
        CoordTransformTuple const& coord_transforms,
        double const TOL = 1e-14)
{
    (check_spline<Patches>(
             random_values_x,
             random_values_y,
             splines_l.template get<Patches>(),
             spline_g,
             std::get<ddc::type_seq_rank_v<Patches, PatchOrdering>>(evaluators_l),
             evaluator_g,
             std::get<ddc::type_seq_rank_v<Patches, PatchOrdering>>(coord_transforms),
             TOL),
     ...);
};

} // end namespace



TEST_F(MultipatchSplineBuilder2DTest, CheckForPeriodicAndGrevilleBC)
{
    IdxRange<GridX<1>, GridY<1>> const& idx_range_xy1 = this->idx_range_xy1;
    IdxRange<GridX<2>, GridY<2>> const& idx_range_xy2 = this->idx_range_xy2;
    IdxRange<GridX<3>, GridY<3>> const& idx_range_xy3 = this->idx_range_xy3;
    IdxRange<GridX<4>, GridY<4>> const& idx_range_xy4 = this->idx_range_xy4;
    IdxRange<GridX<5>, GridY<5>> const& idx_range_xy5 = this->idx_range_xy5;
    IdxRange<GridX<6>, GridY<6>> const& idx_range_xy6 = this->idx_range_xy6;
    IdxRange<GridX<7>, GridY<7>> const& idx_range_xy7 = this->idx_range_xy7;
    IdxRange<GridX<8>, GridY<8>> const& idx_range_xy8 = this->idx_range_xy8;
    IdxRange<GridX<9>, GridY<9>> const& idx_range_xy9 = this->idx_range_xy9;

    IdxRange<GridXg, GridYg> const& idx_range_xy_g = this->idx_range_xy_g;

    // Instantiate the derivatives calculators ---------------------------------------------------
    constexpr std::size_t nb_chosen_cells = 9;

    // InterfaceDerivCoeffs for interfaces along y (periodic).
    InterfaceDerivCoeffs<Interface_1_2> const
            deriv_coeffs_1_2(idx_range_xy1, idx_range_xy2, nb_chosen_cells);
    InterfaceDerivCoeffs<Interface_2_3> const
            deriv_coeffs_2_3(idx_range_xy2, idx_range_xy3, nb_chosen_cells);
    InterfaceDerivCoeffs<Interface_3_1> const
            deriv_coeffs_3_1(idx_range_xy3, idx_range_xy1, nb_chosen_cells);

    InterfaceDerivCoeffs<Interface_4_5> const
            deriv_coeffs_4_5(idx_range_xy4, idx_range_xy5, nb_chosen_cells);
    InterfaceDerivCoeffs<Interface_5_6> const
            deriv_coeffs_5_6(idx_range_xy5, idx_range_xy6, nb_chosen_cells);
    InterfaceDerivCoeffs<Interface_6_4> const
            deriv_coeffs_6_4(idx_range_xy6, idx_range_xy4, nb_chosen_cells);

    InterfaceDerivCoeffs<Interface_7_8> const
            deriv_coeffs_7_8(idx_range_xy7, idx_range_xy8, nb_chosen_cells);
    InterfaceDerivCoeffs<Interface_8_9> const
            deriv_coeffs_8_9(idx_range_xy8, idx_range_xy9, nb_chosen_cells);
    InterfaceDerivCoeffs<Interface_9_7> const
            deriv_coeffs_9_7(idx_range_xy9, idx_range_xy7, nb_chosen_cells);

    // InterfaceDerivCoeffs for interfaces along x.
    InterfaceDerivCoeffs<Interface_1_4> const
            deriv_coeffs_1_4(idx_range_xy1, idx_range_xy4, nb_chosen_cells);
    InterfaceDerivCoeffs<Interface_4_7> const
            deriv_coeffs_4_7(idx_range_xy4, idx_range_xy7, nb_chosen_cells);

    InterfaceDerivCoeffs<Interface_2_5> const
            deriv_coeffs_2_5(idx_range_xy2, idx_range_xy5, nb_chosen_cells);
    InterfaceDerivCoeffs<Interface_5_8> const
            deriv_coeffs_5_8(idx_range_xy5, idx_range_xy8, nb_chosen_cells);

    InterfaceDerivCoeffs<Interface_3_6> const
            deriv_coeffs_3_6(idx_range_xy3, idx_range_xy6, nb_chosen_cells);
    InterfaceDerivCoeffs<Interface_6_9> const
            deriv_coeffs_6_9(idx_range_xy6, idx_range_xy9, nb_chosen_cells);

    // Instantiate the derivative calculators ----------------------------------------------------
    InterfacesDerivativeCalculator<Connectivity, GridX<1>> deriv_calculator_123(
            this->idx_ranges,
            deriv_coeffs_1_2,
            deriv_coeffs_2_3,
            deriv_coeffs_3_1);
    InterfacesDerivativeCalculator<Connectivity, GridX<4>> deriv_calculator_456(
            this->idx_ranges,
            deriv_coeffs_4_5,
            deriv_coeffs_5_6,
            deriv_coeffs_6_4);
    InterfacesDerivativeCalculator<Connectivity, GridX<7>> deriv_calculator_789(
            this->idx_ranges,
            deriv_coeffs_7_8,
            deriv_coeffs_8_9,
            deriv_coeffs_9_7);
    InterfacesDerivativeCalculator<Connectivity, GridY<1>>
            deriv_calculator_147(this->idx_ranges, deriv_coeffs_1_4, deriv_coeffs_4_7);
    InterfacesDerivativeCalculator<Connectivity, GridY<2>>
            deriv_calculator_258(this->idx_ranges, deriv_coeffs_2_5, deriv_coeffs_5_8);
    InterfacesDerivativeCalculator<Connectivity, GridY<3>>
            deriv_calculator_369(this->idx_ranges, deriv_coeffs_3_6, deriv_coeffs_6_9);

    // Instantiate DerivField ====================================================================
    MPDerivFieldMem functions_and_derivs_alloc(this->idx_ranges);
    MPDerivField functions_and_derivs(functions_and_derivs_alloc);

    // Instantiate the field of global function values. No derivatives needed.
    host_t<DFieldMem<IdxRange<GridXg, GridYg>>> function_g_alloc(idx_range_xy_g);
    host_t<DField<IdxRange<GridXg, GridYg>>> function_g = get_field(function_g_alloc);

    // Initialise the data =======================================================================
    // --- the function values.
    initialise_all_functions<Xg, Yg>(functions_and_derivs, this->coord_transforms);
    initialise_2D_function(function_g);

    constexpr ddc::SplineBuilderClosure SBC_H = ddc::SplineBuilderClosure::HERMITE;
    constexpr ddc::SplineBuilderClosure SBC_G = ddc::SplineBuilderClosure::GREVILLE;

    LocalSplineBuilder<Patch1, SBC_H, SBC_H, SBC_H, SBC_G> builder1(idx_range_xy1);
    LocalSplineBuilder<Patch2, SBC_H, SBC_H, SBC_H, SBC_G> builder2(idx_range_xy2);
    LocalSplineBuilder<Patch3, SBC_H, SBC_H, SBC_H, SBC_G> builder3(idx_range_xy3);
    LocalSplineBuilder<Patch4, SBC_H, SBC_H, SBC_H, SBC_H> builder4(idx_range_xy4);
    LocalSplineBuilder<Patch5, SBC_H, SBC_H, SBC_H, SBC_H> builder5(idx_range_xy5);
    LocalSplineBuilder<Patch6, SBC_H, SBC_H, SBC_H, SBC_H> builder6(idx_range_xy6);
    LocalSplineBuilder<Patch7, SBC_H, SBC_H, SBC_G, SBC_H> builder7(idx_range_xy7);
    LocalSplineBuilder<Patch8, SBC_H, SBC_H, SBC_G, SBC_H> builder8(idx_range_xy8);
    LocalSplineBuilder<Patch9, SBC_H, SBC_H, SBC_G, SBC_H> builder9(idx_range_xy9);

    ddc::NullExtrapolationRule null_bc;
    LocalSplineEvaluator<Patch1> evaluator1(null_bc, null_bc, null_bc, null_bc);
    LocalSplineEvaluator<Patch2> evaluator2(null_bc, null_bc, null_bc, null_bc);
    LocalSplineEvaluator<Patch3> evaluator3(null_bc, null_bc, null_bc, null_bc);
    LocalSplineEvaluator<Patch4> evaluator4(null_bc, null_bc, null_bc, null_bc);
    LocalSplineEvaluator<Patch5> evaluator5(null_bc, null_bc, null_bc, null_bc);
    LocalSplineEvaluator<Patch6> evaluator6(null_bc, null_bc, null_bc, null_bc);
    LocalSplineEvaluator<Patch7> evaluator7(null_bc, null_bc, null_bc, null_bc);
    LocalSplineEvaluator<Patch8> evaluator8(null_bc, null_bc, null_bc, null_bc);
    LocalSplineEvaluator<Patch9> evaluator9(null_bc, null_bc, null_bc, null_bc);

    std::tuple evaluators = std::
            tie(evaluator1,
                evaluator2,
                evaluator3,
                evaluator4,
                evaluator5,
                evaluator6,
                evaluator7,
                evaluator8,
                evaluator9);

    using DerivsToCompute = DerivList<
            DerivPerp<GridX<1>>,
            DerivPerp<GridX<4>>,
            DerivPerp<GridX<7>>,
            DerivPerp<GridY<1>>,
            DerivPerp<GridY<2>>,
            DerivPerp<GridY<3>>,
            CrossDeriv<GridY<1>>,
            CrossDeriv<GridY<2>>,
            CrossDeriv<GridY<3>>>;

    MultipatchSplineBuilder2D<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines1OnPatch,
            BSplines2OnPatch,
            Grid1OnPatch,
            Grid2OnPatch,
            ddc::SplineBuilderClosure::GREVILLE,
            ddc::SplineBuilderClosure::GREVILLE,
            ddc::SplineBuilderClosure::GREVILLE,
            ddc::SplineBuilderClosure::GREVILLE,
            ddc::SplineBuilderClosure::HERMITE,
            Connectivity,
            ddc::SplineSolver::LAPACK,
            DConstFieldOnPatch,
            DerivsToCompute,
            DerivativePrecision::APPROXIMATION,
            //     ddc::TypeSeq<
            //             InterfacesDerivativeCalculator<Connectivity, GridX<1>>,
            //             InterfacesDerivativeCalculator<Connectivity, GridX<4>>,
            //             InterfacesDerivativeCalculator<Connectivity, GridX<7>>,
            //             InterfacesDerivativeCalculator<Connectivity, GridY<1>>,
            //             InterfacesDerivativeCalculator<Connectivity, GridY<2>>,
            //             InterfacesDerivativeCalculator<Connectivity, GridY<3>>>
            >
            builders(
                    builder1,
                    builder2,
                    builder3,
                    builder4,
                    builder5,
                    builder6,
                    builder7,
                    builder8,
                    builder9,
                    deriv_calculator_123,
                    deriv_calculator_456,
                    deriv_calculator_789,
                    deriv_calculator_147,
                    deriv_calculator_258,
                    deriv_calculator_369);


    MultipatchFieldMem<
            SplineCoeffMemOnPatch_2D,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            splines_alloc(builders.batched_basis_idx_ranges(this->idx_ranges));
    MultipatchField<
            SplineCoeffOnPatch_2D,
            Patch1,
            Patch2,
            Patch3,
            Patch4,
            Patch5,
            Patch6,
            Patch7,
            Patch8,
            Patch9>
            splines(splines_alloc);

    builders(splines, functions_and_derivs);

    // Test the values of the derivatives ========================================================
    // --- Define an equivalent global spline.
    // Build global spline representation ---
    host_t<DFieldMem<IdxRange<BSplinesXg, BSplinesYg>>> spline_g_alloc(
            this->builder_g.batched_spline_domain(idx_range_xy_g));
    host_t<DField<IdxRange<BSplinesXg, BSplinesYg>>> spline_g = get_field(spline_g_alloc);
    this->builder_g(spline_g, get_const_field(function_g));

    // Check
    std::vector<double> random_values_x = {0, 1. / 3, 0.5, M_PI / 4, 1};
    std::vector<double> random_values_y = {0, 1. / 3, 0.5, M_PI / 4, 1};
    check_all_splines(
            random_values_x,
            random_values_y,
            splines,
            spline_g,
            evaluators,
            this->evaluator_g,
            this->coord_transforms,
            1e-6);
}
