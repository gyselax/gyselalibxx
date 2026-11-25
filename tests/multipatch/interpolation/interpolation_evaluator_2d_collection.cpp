// SPDX-License-Identifier: MIT
#include <array>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_non_uniform.hpp"
#include "constant_extrapolation_rules_onion.hpp"
#include "i_interpolation_builder.hpp"
#include "i_interpolation_evaluator.hpp"
#include "interpolation_evaluator_2d_collection.hpp"
#include "l_norm_tools.hpp"
#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_evaluator.hpp"
#include "mesh_builder.hpp"
#include "multipatch_field.hpp"
#include "multipatch_field_mem.hpp"
#include "multipatch_lagrange_interpolation.hpp"
#include "multipatch_spline_interpolation.hpp"
#include "multipatch_type.hpp"
#include "nd_identity_interpolation_builder.hpp"
#include "nd_lagrange_evaluator.hpp"
#include "null_extrapolation_rules.hpp"
#include "onion_patch_locator.hpp"
#include "physical_geometry.hpp"
#include "spline_testing_tools.hpp"
#include "types.hpp"

template <class Patch>
struct LagBasisR : NonUniformLagrangeBasis<typename Patch::Dim1, 3>
{
};
template <class Patch>
struct LagBasisTheta : NonUniformLagrangeBasis<typename Patch::Dim2, 3>
{
};


namespace {
template <typename ExecSpace>
using MultipatchSplineRThetaInterpolation = MultipatchSplineInterpolator2D<
        PatchLocator<ExecSpace>,
        ConstantExtrapolationRuleOnion<PatchLocator<ExecSpace>>,
        ExtrapolationRule::Constant_Constant,
        ExtrapolationRule::Periodic,
        SplineBoundaryClosure::Greville_Greville,
        SplineBoundaryClosure::Periodic>;

template <typename ExecSpace>
using MultipatchSplineRThetaEvaluator =
        typename MultipatchSplineRThetaInterpolation<ExecSpace>::EvaluatorType;

class InterpolationEvaluator2DCollectionTest : public MultipatchSplineOnionShapeTest
{
protected:
    // Index ranges
    Patch1::IdxRange1 const idx_range_r1;
    Patch1::IdxRange2 const idx_range_theta1;

    Patch2::IdxRange1 const idx_range_r2;
    Patch2::IdxRange2 const idx_range_theta2;

    Patch1::IdxRange1 const idx_range_eval_r1;
    Patch1::IdxRange12 const idx_range_eval_rtheta1;

    Patch2::IdxRange1 const idx_range_eval_r2;
    Patch2::IdxRange12 const idx_range_eval_rtheta2;

    // Spline representations
    DField<IdxRange<BSplinesR<1>, BSplinesTheta<1>>> const function_1_coef;
    DField<IdxRange<BSplinesR<2>, BSplinesTheta<2>>> const function_2_coef;

    LogicalToPhysicalMapping const to_physical_mapping;
    PhysicalToLogicalMapping const to_logical_mapping;

    // Interpolations on every patch on device and host
    MultipatchSplineRThetaInterpolation<DeviceExecSpace> const interpolation;
    MultipatchSplineRThetaInterpolation<HostExecSpace> const interpolation_host;

    // Evaluators on every patch
    MultipatchSplineRThetaEvaluator<DeviceExecSpace> const& evaluators;
    MultipatchSplineRThetaEvaluator<HostExecSpace> const& evaluators_host;

    // Evaluators on single patch
    SplineRThetaEvaluator<1, DeviceExecSpace> const& evaluator_1;
    SplineRThetaEvaluator<2, DeviceExecSpace> const& evaluator_2;

    SplineRThetaEvaluator<1, HostExecSpace> const& evaluator_1_host;
    SplineRThetaEvaluator<2, HostExecSpace> const& evaluator_2_host;

public:
    InterpolationEvaluator2DCollectionTest()
        // Idx ranges
        : idx_range_r1(SplineInterpPointsR<1>::get_domain<GridR<1>>())
        , idx_range_theta1(SplineInterpPointsTheta<1>::get_domain<GridTheta<1>>())
        , idx_range_r2(SplineInterpPointsR<2>::get_domain<GridR<2>>())
        , idx_range_theta2(SplineInterpPointsTheta<2>::get_domain<GridTheta<2>>())
        , idx_range_eval_r1(idx_range_r1.remove_last(Patch1::IdxStep1(1)))
        , idx_range_eval_rtheta1(
                  idx_range_eval_r1,
                  idx_range_theta1.remove_last(Patch1::IdxStep2(1)))
        , idx_range_eval_r2(idx_range_r2.remove_last(Patch2::IdxStep1(1)))
        , idx_range_eval_rtheta2(
                  idx_range_eval_r2,
                  idx_range_theta2.remove_last(Patch2::IdxStep2(1)))
        // Local splines on device and host
        , function_1_coef(get_field(function_1_coef_alloc))
        , function_2_coef(get_field(function_2_coef_alloc))
        , to_physical_mapping()
        , to_logical_mapping()
        , interpolation(
                  all_idx_ranges,
                  PatchLocator<
                          DeviceExecSpace>(all_idx_ranges, to_physical_mapping, to_logical_mapping),
                  ConstantExtrapolationRuleOnion<PatchLocator<DeviceExecSpace>>(r1_min, r2_max))
        , interpolation_host(
                  all_idx_ranges,
                  PatchLocator<
                          HostExecSpace>(all_idx_ranges, to_physical_mapping, to_logical_mapping),
                  ConstantExtrapolationRuleOnion<PatchLocator<HostExecSpace>>(r1_min, r2_max))
        , evaluators(interpolation.get_evaluator())
        , evaluators_host(interpolation_host.get_evaluator())
        , evaluator_1(interpolation.get_interpolation<Patch1>().get_evaluator())
        , evaluator_2(interpolation.get_interpolation<Patch2>().get_evaluator())
        , evaluator_1_host(interpolation_host.get_interpolation<Patch1>().get_evaluator())
        , evaluator_2_host(interpolation_host.get_interpolation<Patch2>().get_evaluator())
    {
    }
};

template <class GridR, class GridTheta>
void set_eval_points_2D(
        Field<Coord<typename GridR::continuous_dimension_type,
                    typename GridTheta::continuous_dimension_type>,
              IdxRange<GridR, GridTheta>> eval_points)
{
    using CoordRTheta
            = Coord<typename GridR::continuous_dimension_type,
                    typename GridTheta::continuous_dimension_type>;
    using CoordR = Coord<typename GridR::continuous_dimension_type>;
    using CoordTheta = Coord<typename GridTheta::continuous_dimension_type>;

    IdxRange<GridR, GridTheta> const idx_range_eval = get_idx_range(eval_points);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_eval,
            KOKKOS_LAMBDA(Idx<GridR, GridTheta> const idx) {
                CoordRTheta const coord = ddc::coordinate(idx);
                CoordR const delta_r = distance_at_right(Idx<GridR>(idx));
                CoordTheta const delta_theta = distance_at_right(Idx<GridTheta>(idx));
                eval_points(idx) = coord + CoordRTheta(delta_r / 2., delta_theta / 2.);
            });
}

/// Initialise f(r, theta) = r^power sin(theta) on the interpolation grid of a patch.
template <class Patch>
void initialise_function(DField<typename Patch::IdxRange12> values, int const power)
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(values),
            KOKKOS_LAMBDA(typename Patch::Idx12 const idx) {
                double const r = ddc::coordinate(typename Patch::Idx1(idx));
                double const theta = ddc::coordinate(typename Patch::Idx2(idx));
                values(idx) = Kokkos::pow(r, power) * Kokkos::sin(theta);
            });
}

// Checking functions ----------------------------------------------------------------------------
// One function per method of InterpolationEvaluator2DCollection tested.
// Test operator() ...............................................................................
void test_operator_assignment(
        MultipatchSplineRThetaEvaluator<DeviceExecSpace> const& evaluators,
        SplineRThetaEvaluator<1, DeviceExecSpace> const& single_evaluator_1,
        SplineRThetaEvaluator<2, DeviceExecSpace> const& single_evaluator_2,
        Kokkos::View<Coord<R, Theta>*, Kokkos::DefaultExecutionSpace::memory_space> const&
                eval_coords,
        MultipatchField<ConstSplineCoeffOnPatch_2D, Patch1, Patch2> const& splines,
        ConstSplineCoeffOnPatch_2D<Patch1> const& spline_patch_1,
        ConstSplineCoeffOnPatch_2D<Patch2> const& spline_patch_2)
{
    double max_error = 0;
    Kokkos::parallel_reduce(
            eval_coords.extent(0),
            KOKKOS_LAMBDA(int i, double& err) {
                double const eval_function = evaluators(eval_coords(i), splines);
                double expected_function;
                if (Coord<R>(eval_coords(i)) < ddc::discrete_space<BSplinesR<1>>().rmax()) {
                    expected_function = single_evaluator_1(eval_coords(i), spline_patch_1);
                } else {
                    expected_function = single_evaluator_2(eval_coords(i), spline_patch_2);
                }
                err = Kokkos::max(Kokkos::abs(eval_function - expected_function), err);
                err = 0;
            },
            Kokkos::Max<double>(max_error));
    EXPECT_LE(max_error, 1e-12);
}

// Test deriv() on device and on host ............................................................
template <class... DerivDims, class Evaluator, class IdxRangeBS, class MemorySpace>
KOKKOS_FUNCTION double call_correct_operator(
        Idx<DerivDims...> const& deriv_order,
        Evaluator const& evaluator,
        Coord<R, Theta> const& eval_coord,
        DConstField<IdxRangeBS, MemorySpace> const& spline)
{
    static_assert(
            (ddc::in_tags_v<DerivDims, ddc::TypeSeq<ddc::Deriv<R>, ddc::Deriv<Theta>>> && ...),
            "Please provide a deriv_order indexed on the derivative dimensions of the patches "
            "(same dimensions for all the patches).");
    using IdxDeriv = Idx<DerivDims...>;
    const Idx<ddc::Deriv<R>> dr_order0(0);
    const Idx<ddc::Deriv<Theta>> dtheta_order0(0);
    if constexpr (std::is_same_v<IdxDeriv, Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>>) {
        Idx<ddc::Deriv<R>> dr_order(deriv_order);
        Idx<ddc::Deriv<Theta>> dtheta_order(deriv_order);
        if ((dr_order == dr_order0) && (dtheta_order == dtheta_order0)) {
            return evaluator(eval_coord, spline);
        } else if (dr_order == dr_order0) {
            return evaluator.deriv(dtheta_order, eval_coord, spline);
        } else if (dtheta_order == dtheta_order0) {
            return evaluator.deriv(dr_order, eval_coord, spline);
        } else {
            return evaluator.deriv(deriv_order, eval_coord, spline);
        }
    } else if constexpr (std::is_same_v<IdxDeriv, Idx<ddc::Deriv<R>>>) {
        if (deriv_order == dr_order0) {
            return evaluator(eval_coord, spline);
        } else {
            return evaluator.deriv(deriv_order, eval_coord, spline);
        }
    } else {
        static_assert(std::is_same_v<IdxDeriv, Idx<ddc::Deriv<Theta>>>);
        if (deriv_order == dtheta_order0) {
            return evaluator(eval_coord, spline);
        } else {
            return evaluator.deriv(deriv_order, eval_coord, spline);
        }
    }
}

template <template <typename P> class SplineTypeOnPatch, class IdxDeriv, class ExecSpace>
void test_deriv_at_order(
        IdxDeriv const& deriv_order,
        MultipatchSplineRThetaEvaluator<ExecSpace> const& evaluators,
        SplineRThetaEvaluator<1, ExecSpace> const& single_evaluator_1,
        SplineRThetaEvaluator<2, ExecSpace> const& single_evaluator_2,
        Kokkos::View<Coord<R, Theta>*, typename ExecSpace::memory_space> const& eval_coords,
        MultipatchField<SplineTypeOnPatch, Patch1, Patch2> const& splines,
        SplineTypeOnPatch<Patch1> spline_patch_1,
        SplineTypeOnPatch<Patch2> spline_patch_2,
        double const TOL)
{
    double max_error = 0;
    Kokkos::parallel_reduce(
            Kokkos::RangePolicy<ExecSpace>(0, eval_coords.extent(0)),
            KOKKOS_LAMBDA(int i, double& err) {
                double const eval_function = evaluators.deriv(deriv_order, eval_coords(i), splines);
                double expected_function;
                if (Coord<R>(eval_coords(i)) < ddc::discrete_space<BSplinesR<1>>().rmax()) {
                    expected_function = call_correct_operator(
                            deriv_order,
                            single_evaluator_1,
                            eval_coords(i),
                            spline_patch_1);
                } else {
                    expected_function = call_correct_operator(
                            deriv_order,
                            single_evaluator_2,
                            eval_coords(i),
                            spline_patch_2);
                }
                err = Kokkos::max(Kokkos::abs(eval_function - expected_function), err);
            },
            Kokkos::Max<double>(max_error));
    EXPECT_LE(max_error, TOL);
};

template <template <typename P> class SplineTypeOnPatch, class ExecSpace>
void test_deriv(
        MultipatchSplineRThetaEvaluator<ExecSpace> const& evaluators,
        SplineRThetaEvaluator<1, ExecSpace> const& single_evaluator_1,
        SplineRThetaEvaluator<2, ExecSpace> const& single_evaluator_2,
        Kokkos::View<Coord<R, Theta>*, typename ExecSpace::memory_space> const& eval_coords,
        MultipatchField<SplineTypeOnPatch, Patch1, Patch2> const& splines,
        SplineTypeOnPatch<Patch1> const& spline_patch_1,
        SplineTypeOnPatch<Patch2> const& spline_patch_2)
{
    // Test deriv on dimension R.
    for (int i_dr(0); i_dr <= 3; i_dr++) {
        Idx<ddc::Deriv<R>> deriv_order(i_dr);
        test_deriv_at_order<SplineTypeOnPatch>(
                deriv_order,
                evaluators,
                single_evaluator_1,
                single_evaluator_2,
                eval_coords,
                splines,
                spline_patch_1,
                spline_patch_2,
                1e-11);
    }
    // Test deriv on dimension Theta.
    for (int i_dtheta(0); i_dtheta <= 3; i_dtheta++) {
        Idx<ddc::Deriv<Theta>> deriv_order(i_dtheta);
        test_deriv_at_order<SplineTypeOnPatch>(
                deriv_order,
                evaluators,
                single_evaluator_1,
                single_evaluator_2,
                eval_coords,
                splines,
                spline_patch_1,
                spline_patch_2,
                1e-11);
    }
    // Test deriv on dimensions (R, Theta).
    for (int i_dr(0); i_dr <= 3; i_dr++) {
        for (int i_dtheta(0); i_dtheta <= 3; i_dtheta++) {
            Idx<ddc::Deriv<R>, ddc::Deriv<Theta>> deriv_order(i_dr, i_dtheta);
            test_deriv_at_order<SplineTypeOnPatch>(
                    deriv_order,
                    evaluators,
                    single_evaluator_1,
                    single_evaluator_2,
                    eval_coords,
                    splines,
                    spline_patch_1,
                    spline_patch_2,
                    5e-10);
        }
    }
};

template <class PatchP>
void check_near(
        DConstFieldOnPatch_host<PatchP> expected_host,
        DConstFieldOnPatch_host<PatchP> eval_host,
        double const TOL)
{
    ddc::host_for_each(get_idx_range(expected_host), [&](typename PatchP::Idx12 const idx) {
        EXPECT_NEAR(expected_host(idx), eval_host(idx), TOL);
    });
};

template <typename PatchP, typename GridDeriv>
void check_deriv_match(
        Idx<ddc::Deriv<typename GridDeriv::continuous_dimension_type>> deriv_order,
        typename PatchP::IdxRange12 const& idx_range_eval_rtheta,
        DerivFieldOnPatch<PatchP> eval_function_and_derivs,
        SplineRThetaEvaluator<
                ddc::type_seq_rank_v<PatchP, ddc::TypeSeq<Patch1, Patch2>> + 1,
                DeviceExecSpace> const& evaluator,
        DConstField<typename PatchP::IdxRangeBS12> spline)
{
    using GridEdge = std::conditional_t<
            std::is_same_v<GridDeriv, typename PatchP::Grid1>,
            typename PatchP::Grid2,
            typename PatchP::Grid1>;

    IdxRange<GridDeriv> idx_range_eval_1d(idx_range_eval_rtheta);
    IdxRange<GridEdge> idx_range_eval_edge(idx_range_eval_rtheta);

    for (Idx<GridDeriv> idx_1d :
         std::vector<Idx<GridDeriv>> {idx_range_eval_1d.front(), idx_range_eval_1d.back()}) {
        IdxRange<GridDeriv> const idx_range_1d(idx_1d, IdxStep<GridDeriv>(1));
        typename PatchP::IdxRange12 idx_range_deriv(idx_range_1d, idx_range_eval_edge);

        // Evaluation points on the edge.
        CoordFieldMemOnPatch<PatchP> eval_points_deriv_alloc(idx_range_deriv);
        CoordFieldOnPatch<PatchP> eval_points_deriv(eval_points_deriv_alloc);

        set_eval_points_2D<typename PatchP::Grid1, typename PatchP::Grid2>(eval_points_deriv);

        // Expected functions.
        DFieldMemOnPatch<PatchP> expected_deriv_alloc(idx_range_deriv);
        DFieldOnPatch<PatchP> expected_deriv(expected_deriv_alloc);

        if (deriv_order == Idx<ddc::Deriv<typename GridDeriv::continuous_dimension_type>>(0)) {
            evaluator(expected_deriv, get_const_field(eval_points_deriv), spline);
        } else {
            evaluator
                    .deriv(deriv_order, expected_deriv, get_const_field(eval_points_deriv), spline);
        }

        // Compare the evaluated derivatives with the expected derivatives.
        // --- define DField for ddc::create_mirror_and_copy().
        DFieldMemOnPatch<PatchP> eval_deriv_alloc(idx_range_deriv);
        DFieldOnPatch<PatchP> eval_deriv(eval_deriv_alloc);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(eval_deriv),
                KOKKOS_LAMBDA(typename PatchP::Idx12 const idx) {
                    eval_deriv(idx) = eval_function_and_derivs(deriv_order, idx);
                });

        // --- put the functions on host
        auto const expected_deriv_host = ddc::create_mirror_and_copy(expected_deriv);
        auto const eval_deriv_host = ddc::create_mirror_and_copy(eval_deriv);

        // --- check errors.
        check_near<PatchP>(get_field(expected_deriv_host), get_field(eval_deriv_host), 1e-11);
    }
};

template <typename PatchP>
void check_cross_deriv_match(
        Idx<ddc::Deriv<typename PatchP::Dim1>, ddc::Deriv<typename PatchP::Dim2>> deriv_order,
        typename PatchP::IdxRange12 const& idx_range_eval_rtheta,
        DerivFieldOnPatch<PatchP> eval_function_and_derivs,
        SplineRThetaEvaluator<
                ddc::type_seq_rank_v<PatchP, ddc::TypeSeq<Patch1, Patch2>> + 1,
                DeviceExecSpace> const& evaluator,
        DConstField<typename PatchP::IdxRangeBS12> spline)
{
    typename PatchP::IdxRange1 idx_range_eval_1(idx_range_eval_rtheta);
    typename PatchP::IdxRange2 idx_range_eval_2(idx_range_eval_rtheta);

    for (typename PatchP::Idx1 idx_1 :
         std::vector<typename PatchP::Idx1> {idx_range_eval_1.front(), idx_range_eval_1.back()}) {
        for (typename PatchP::Idx2 idx_2 : std::vector<
                     typename PatchP::Idx2> {idx_range_eval_2.front(), idx_range_eval_2.back()}) {
            typename PatchP::IdxRange1 const idx_range_1(idx_1, typename PatchP::IdxStep1(1));
            typename PatchP::IdxRange2 const idx_range_2(idx_2, typename PatchP::IdxStep2(1));
            typename PatchP::IdxRange12 idx_range_deriv(idx_range_1, idx_range_2);

            // Evaluation points on the edge.
            CoordFieldMemOnPatch<PatchP> eval_points_cross_deriv_alloc(idx_range_deriv);
            CoordFieldOnPatch<PatchP> eval_points_cross_deriv(eval_points_cross_deriv_alloc);

            set_eval_points_2D<typename PatchP::Grid1, typename PatchP::Grid2>(
                    eval_points_cross_deriv);

            // Expected functions.
            DFieldMemOnPatch<PatchP> expected_cross_deriv_alloc(idx_range_deriv);
            DFieldOnPatch<PatchP> expected_cross_deriv(expected_cross_deriv_alloc);

            evaluator
                    .deriv(deriv_order,
                           expected_cross_deriv,
                           get_const_field(eval_points_cross_deriv),
                           spline);

            // Compare the evaluated derivatives with the expected derivatives.
            // --- define DField for ddc::create_mirror_and_copy().
            DFieldMemOnPatch<PatchP> eval_cross_deriv_alloc(idx_range_deriv);
            DFieldOnPatch<PatchP> eval_cross_deriv(eval_cross_deriv_alloc);

            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    get_idx_range(eval_cross_deriv),
                    KOKKOS_LAMBDA(typename PatchP::Idx12 const idx) {
                        eval_cross_deriv(idx) = eval_function_and_derivs(deriv_order, idx);
                    });

            // --- put the functions on host
            auto const expected_cross_deriv_host
                    = ddc::create_mirror_and_copy(expected_cross_deriv);
            auto const eval_cross_deriv_host = ddc::create_mirror_and_copy(eval_cross_deriv);

            // --- check errors.
            check_near<PatchP>(
                    get_field(expected_cross_deriv_host),
                    get_field(eval_cross_deriv_host),
                    1e-10);
        }
    }
};

// Non-spline (Lagrange) evaluators ------------------------------------------------------------
template <class Patch>
using LagrangeEvaluatorR = LagrangeEvaluator<
        DeviceExecSpace,
        typename DeviceExecSpace::memory_space,
        double,
        LagBasisR<Patch>,
        typename Patch::Grid1,
        ddc::NullExtrapolationRule,
        ddc::NullExtrapolationRule>;

template <class Patch>
using LagrangeEvaluatorTheta = LagrangeEvaluator<
        DeviceExecSpace,
        typename DeviceExecSpace::memory_space,
        double,
        LagBasisTheta<Patch>,
        typename Patch::Grid2,
        ddc::PeriodicExtrapolationRule<typename Patch::Dim2>,
        ddc::PeriodicExtrapolationRule<typename Patch::Dim2>>;

template <class Patch>
using LagrangeEvaluatorRTheta
        = NDLagrangeEvaluator<LagrangeEvaluatorR<Patch>, LagrangeEvaluatorTheta<Patch>>;

template <class Patch>
using LagrangeBuilderRTheta = NDIdentityInterpolationBuilder<
        DeviceExecSpace,
        typename DeviceExecSpace::memory_space,
        double,
        typename Patch::IdxRange12,
        IdxRange<LagBasisR<Patch>, LagBasisTheta<Patch>>>;

template <class Patch>
using IdxRangeLagrangeCoeff =
        typename InterpolationEvaluatorTraits<LagrangeEvaluatorRTheta<Patch>>::coeff_idx_range_type;

template <class Patch>
using LagrangeCoeffMemOnPatch = DFieldMem<IdxRangeLagrangeCoeff<Patch>>;

template <class Patch>
using LagrangeConstCoeffOnPatch = DConstField<IdxRangeLagrangeCoeff<Patch>>;

using MultipatchLagrangeInterpolation = MultipatchLagrangeInterpolator2D<
        PatchLocator<DeviceExecSpace>,
        NullExtrapolationRule,
        LagBasisR,
        LagBasisTheta,
        ExtrapolationRule::Null_Null,
        ExtrapolationRule::Periodic>;

using LagrangeEvaluatorCollection = MultipatchLagrangeInterpolation::EvaluatorType;

/**
 * Initialise the (non-uniform) grids of a patch with uniformly spaced points and the
 * Lagrange bases defined on these grids. The second dimension is periodic on [0, 2pi].
 * @return The index range of the interpolation points (without the last periodic point).
 */
template <class Patch>
typename Patch::IdxRange12 init_lagrange_patch(
        typename Patch::Coord1 r_min,
        typename Patch::Coord1 r_max,
        std::size_t r_ncells,
        std::size_t theta_ncells)
{
    using Grid1 = typename Patch::Grid1;
    using Grid2 = typename Patch::Grid2;
    typename Patch::Coord2 const theta_min(0.0);
    typename Patch::Coord2 const theta_max(2 * M_PI);
    ddc::init_discrete_space<Grid1>(
            build_uniform_break_points(r_min, r_max, IdxStep<Grid1>(r_ncells)));
    ddc::init_discrete_space<Grid2>(
            build_uniform_break_points(theta_min, theta_max, IdxStep<Grid2>(theta_ncells)));
    IdxRange<Grid1> idx_range1(Idx<Grid1>(0), IdxStep<Grid1>(r_ncells + 1));
    IdxRange<Grid2> idx_range2(Idx<Grid2>(0), IdxStep<Grid2>(theta_ncells + 1));
    ddc::init_discrete_space<LagBasisR<Patch>>(idx_range1);
    ddc::init_discrete_space<LagBasisTheta<Patch>>(idx_range2);
    // The last point in the periodic dimension is the same as the first point.
    return typename Patch::IdxRange12(idx_range1, idx_range2.remove_last(IdxStep<Grid2>(1)));
}

/**
 * Build the Lagrange coefficients of the function f(r, theta) = r^power sin(theta)
 * on the given patch.
 */
template <class Patch>
LagrangeCoeffMemOnPatch<Patch> build_lagrange_coeffs(
        typename Patch::IdxRange12 const& idx_range,
        int const power)
{
    DFieldMem<typename Patch::IdxRange12> values_alloc(idx_range);
    DField<typename Patch::IdxRange12> values = get_field(values_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(typename Patch::Idx12 const idx) {
                typename Patch::Idx1 const idx_r(idx);
                typename Patch::Idx2 const idx_theta(idx);
                double const r = ddc::coordinate(idx_r);
                double const theta = ddc::coordinate(idx_theta);
                values(idx) = Kokkos::pow(r, power) * Kokkos::sin(theta);
            });

    LagrangeBuilderRTheta<Patch> const builder;
    static_assert(std::is_same_v<
                  typename InterpolationBuilderTraits<LagrangeBuilderRTheta<Patch>>::
                          template batched_basis_idx_range_type<typename Patch::IdxRange12>,
                  IdxRangeLagrangeCoeff<Patch>>);
    LagrangeCoeffMemOnPatch<Patch> coeffs_alloc(batched_basis_idx_range(builder, idx_range));
    builder(get_field(coeffs_alloc), get_const_field(values_alloc));
    return coeffs_alloc;
}

/// Create the Lagrange evaluator of a patch.
template <class Patch>
LagrangeEvaluatorRTheta<Patch> make_lagrange_evaluator()
{
    ddc::NullExtrapolationRule const null_extrap;
    ddc::PeriodicExtrapolationRule<typename Patch::Dim2> const periodic_extrap;
    LagrangeEvaluatorR<Patch> const evaluator_r(null_extrap, null_extrap);
    LagrangeEvaluatorTheta<Patch> const evaluator_theta(periodic_extrap, periodic_extrap);
    return LagrangeEvaluatorRTheta<Patch>(evaluator_r, evaluator_theta);
}

/**
 * @brief Google test class to define the Lagrange interpolations on an Onion shape
 * domain with 2 patches.
 */
class InterpolationEvaluator2DCollectionLagrangeTest : public ::testing::Test
{
protected:
    Patch1::IdxRange12 const idx_range_rtheta1;
    Patch2::IdxRange12 const idx_range_rtheta2;

    // Lagrange coefficients of f(r, theta) = r sin(theta) on Patch1
    // and f(r, theta) = r^2 sin(theta) on Patch2.
    LagrangeCoeffMemOnPatch<Patch1> const function_1_coef_alloc;
    LagrangeCoeffMemOnPatch<Patch2> const function_2_coef_alloc;
    MultipatchField<LagrangeConstCoeffOnPatch, Patch1, Patch2> const coeffs;

    MultipatchIdxRange const all_idx_ranges;
    LogicalToPhysicalMapping const to_physical_mapping;
    PhysicalToLogicalMapping const to_logical_mapping;

    // Interpolations on every patch (the NullExtrapolationRule is default-constructed)
    MultipatchLagrangeInterpolation const interpolation;

    // Evaluators on every patch
    LagrangeEvaluatorCollection const& evaluators;

    // Evaluators on single patch
    LagrangeEvaluatorRTheta<Patch1> const& evaluator_1;
    LagrangeEvaluatorRTheta<Patch2> const& evaluator_2;

public:
    InterpolationEvaluator2DCollectionLagrangeTest()
        : idx_range_rtheta1(
                init_lagrange_patch<Patch1>(Patch1::Coord1(0.25), Patch1::Coord1(1.0), 16, 10))
        , idx_range_rtheta2(
                  init_lagrange_patch<Patch2>(Patch2::Coord1(1.0), Patch2::Coord1(2.0), 8, 12))
        , function_1_coef_alloc(build_lagrange_coeffs<Patch1>(idx_range_rtheta1, 1))
        , function_2_coef_alloc(build_lagrange_coeffs<Patch2>(idx_range_rtheta2, 2))
        , coeffs(get_const_field(function_1_coef_alloc), get_const_field(function_2_coef_alloc))
        , all_idx_ranges(idx_range_rtheta1, idx_range_rtheta2)
        , to_physical_mapping()
        , to_logical_mapping()
        , interpolation(
                  all_idx_ranges,
                  PatchLocator<
                          DeviceExecSpace>(all_idx_ranges, to_physical_mapping, to_logical_mapping))
        , evaluators(interpolation.get_evaluator())
        , evaluator_1(interpolation.get_interpolation<Patch1>().get_evaluator())
        , evaluator_2(interpolation.get_interpolation<Patch2>().get_evaluator())
    {
    }
};

/**
 * Check the collection against the evaluator of each patch on fields of coordinates
 * (each coordinate is located on its storing patch).
 */
template <class... DerivDims>
void check_lagrange_on_coord_field(
        Idx<DerivDims...> const& deriv_order,
        LagrangeEvaluatorCollection const& evaluators,
        LagrangeEvaluatorRTheta<Patch1> const& evaluator_1,
        LagrangeEvaluatorRTheta<Patch2> const& evaluator_2,
        MultipatchField<CoordConstFieldOnPatch, Patch1, Patch2> const& eval_points,
        MultipatchField<LagrangeConstCoeffOnPatch, Patch1, Patch2> const& coeffs)
{
    Patch1::IdxRange12 const idx_range_eval_1 = get_idx_range(eval_points.get<Patch1>());
    Patch2::IdxRange12 const idx_range_eval_2 = get_idx_range(eval_points.get<Patch2>());

    MultipatchType<IdxRangeOnPatch, Patch1, Patch2> idx_range(idx_range_eval_1, idx_range_eval_2);

    MultipatchFieldMem<DFieldMemOnPatch, Patch1, Patch2> eval_alloc(idx_range);
    MultipatchFieldMem<DFieldMemOnPatch, Patch1, Patch2> expected_alloc(idx_range);

    MultipatchField<DFieldOnPatch, Patch1, Patch2> const eval_values = get_field(eval_alloc);

    if constexpr (sizeof...(DerivDims) == 0) {
        evaluators(eval_values, eval_points, coeffs);
        evaluator_1(expected_alloc.get<Patch1>(), eval_points.get<Patch1>(), coeffs.get<Patch1>());
        evaluator_2(expected_alloc.get<Patch2>(), eval_points.get<Patch2>(), coeffs.get<Patch2>());
    } else {
        evaluators.deriv(deriv_order, eval_values, eval_points, coeffs);
        evaluator_1
                .deriv(deriv_order,
                       expected_alloc.get<Patch1>(),
                       eval_points.get<Patch1>(),
                       coeffs.get<Patch1>());
        evaluator_2
                .deriv(deriv_order,
                       expected_alloc.get<Patch2>(),
                       eval_points.get<Patch2>(),
                       coeffs.get<Patch2>());
    }

    EXPECT_LE(
            error_norm_inf(
                    Kokkos::DefaultExecutionSpace(),
                    get_const_field(eval_alloc.get<Patch1>()),
                    get_const_field(expected_alloc.get<Patch1>())),
            1e-12);
    EXPECT_LE(
            error_norm_inf(
                    Kokkos::DefaultExecutionSpace(),
                    get_const_field(eval_alloc.get<Patch2>()),
                    get_const_field(expected_alloc.get<Patch2>())),
            1e-12);
}

/**
 * Check the collection against the evaluator of the patch where each coordinate is
 * physically located. Outside the domain the null extrapolation rule is expected.
 */
void check_lagrange_on_single_coords(
        LagrangeEvaluatorCollection const& evaluators,
        LagrangeEvaluatorRTheta<Patch1> const& evaluator_1,
        LagrangeEvaluatorRTheta<Patch2> const& evaluator_2,
        Kokkos::View<Coord<R, Theta>*, Kokkos::DefaultExecutionSpace::memory_space> const&
                eval_coords,
        MultipatchField<LagrangeConstCoeffOnPatch, Patch1, Patch2> const& coeffs)
{
    LagrangeConstCoeffOnPatch<Patch1> const coeffs_1 = coeffs.get<Patch1>();
    LagrangeConstCoeffOnPatch<Patch2> const coeffs_2 = coeffs.get<Patch2>();
    // The interface between the patches is at r = 1 and the outer boundary at r = 2.
    double max_error = 0;
    Kokkos::parallel_reduce(
            eval_coords.extent(0),
            KOKKOS_LAMBDA(int i, double& err) {
                Coord<R, Theta> const coord = eval_coords(i);
                double const eval_function = evaluators(coord, coeffs);
                double expected_function = 0.0;
                double const r = ddc::get<R>(coord);
                if (r < 1.0) {
                    expected_function = evaluator_1(coord, coeffs_1);
                } else if (r <= 2.0) {
                    expected_function = evaluator_2(coord, coeffs_2);
                }
                err = Kokkos::max(Kokkos::abs(eval_function - expected_function), err);
            },
            Kokkos::Max<double>(max_error));
    EXPECT_LE(max_error, 1e-12);
}

} // namespace


/* -----------------------------------------------------------------------------------------------
    Test operator() for a single coordinate on host.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionTest, HostEvaluateOnSingleCoord)
{
    // Compare the evaluated functions with the expected functions.
    // Evaluation points
    // --- patch 1
    typename Patch1::Coord12 const eval_coord_P1_on_patch_1(0.5, 0.0);
    typename Patch1::Coord12 const eval_coord_P1_on_patch_2(1.1, M_PI);

    // --- patch 2
    typename Patch2::Coord12 const eval_coord_P2_on_patch_1(0.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 const eval_coord_P2_on_patch_2(1.9, 4 * M_PI);
    typename Patch2::Coord12 const eval_coord_P2_outside(2.3, -M_PI);

    auto const function_1_coef_host = ddc::create_mirror_and_copy(function_1_coef);
    auto const function_2_coef_host = ddc::create_mirror_and_copy(function_2_coef);
    MultipatchField<ConstSplineCoeffOnPatch_2D_host, Patch1, Patch2> const
            splines_host(get_field(function_1_coef_host), get_field(function_2_coef_host));

    host_t<DConstField<IdxRange<BSplinesR<1>, BSplinesTheta<1>>>> const const_function_1_coef_host(
            get_const_field(function_1_coef_host));
    host_t<DConstField<IdxRange<BSplinesR<2>, BSplinesTheta<2>>>> const const_function_2_coef_host(
            get_const_field(function_2_coef_host));

    // Check assignment operator on host.
    double eval_function = evaluators_host(eval_coord_P1_on_patch_1, splines_host);
    double expected_function
            = evaluator_1_host(eval_coord_P1_on_patch_1, const_function_1_coef_host);
    EXPECT_NEAR(eval_function, expected_function, 1e-12);

    eval_function = evaluators_host(eval_coord_P1_on_patch_2, splines_host);
    expected_function = evaluator_2_host(eval_coord_P1_on_patch_2, const_function_2_coef_host);
    EXPECT_NEAR(eval_function, expected_function, 1e-12);

    eval_function = evaluators_host(eval_coord_P2_on_patch_1, splines_host);
    expected_function = evaluator_1_host(eval_coord_P2_on_patch_1, const_function_1_coef_host);
    EXPECT_NEAR(eval_function, expected_function, 1e-12);

    eval_function = evaluators_host(eval_coord_P2_on_patch_2, splines_host);
    expected_function = evaluator_2_host(eval_coord_P2_on_patch_2, const_function_2_coef_host);
    EXPECT_NEAR(eval_function, expected_function, 1e-12);

    eval_function = evaluators_host(eval_coord_P2_outside, splines_host);
    expected_function = evaluator_2_host(eval_coord_P2_outside, const_function_2_coef_host);
    EXPECT_NEAR(eval_function, expected_function, 1e-12);

    /*
        REMARK: here Patch1::Coord12 == Patch2::Coord12 == Coord<R, Theta>.
            So we do not need conversion to apply the spline evaluators of the
            patch 1 a coordinate from the patch 2, and vice versa.
    */
}


/* -----------------------------------------------------------------------------------------------
    Test operator() for a single coordinate on device.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionTest, DeviceEvaluateOnSingleCoord)
{
    // Compare the evaluated functions with the expected functions.
    // Evaluation points
    // --- patch 1
    typename Patch1::Coord12 const eval_coord_P1_on_patch_1(0.5, 0.0);
    typename Patch1::Coord12 const eval_coord_P1_on_patch_2(1.1, M_PI);

    // --- patch 2
    typename Patch2::Coord12 const eval_coord_P2_on_patch_1(0.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 const eval_coord_P2_on_patch_2(1.9, 4 * M_PI);
    typename Patch2::Coord12 const eval_coord_P2_outside(2.3, -M_PI);

    // Coordinates and patch indices on device.
    Kokkos::View<Coord<R, Theta>*, Kokkos::DefaultExecutionSpace::memory_space> coords("coords", 5);
    Kokkos::View<Coord<R, Theta>*, Kokkos::HostSpace> coords_host("coords_host", 5);

    coords_host(0) = eval_coord_P1_on_patch_1;
    coords_host(1) = eval_coord_P1_on_patch_2;
    coords_host(2) = eval_coord_P2_on_patch_1;
    coords_host(3) = eval_coord_P2_on_patch_2;
    coords_host(4) = eval_coord_P2_outside;

    Kokkos::deep_copy(coords, coords_host);

    // Check assignment operator on device.
    test_operator_assignment(
            evaluators,
            evaluator_1,
            evaluator_2,
            coords,
            splines,
            get_const_field(function_1_coef),
            get_const_field(function_2_coef));

    /*
        REMARK: here Patch1::Coord12 == Patch2::Coord12 == Coord<R, Theta>.
            So we do not need conversion to apply the spline evaluators of the
            patch 1 a coordinate from the patch 2, and vice versa.
    */
}


/* -----------------------------------------------------------------------------------------------
    Test operator() for fields of coordinates called from host but stored on device.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionTest, EvaluateOnCoordField)
{
    // Evaluation points
    // --- patch 1
    Patch1::IdxRange12 const idx_range_eval_periodic_rtheta1(idx_range_eval_r1, idx_range_theta1);
    FieldMem<Patch1::Coord12, Patch1::IdxRange12> eval_points_1_alloc(
            idx_range_eval_periodic_rtheta1);
    Field<Patch1::Coord12, Patch1::IdxRange12> eval_points_1 = get_field(eval_points_1_alloc);

    // --- patch 2
    Patch2::IdxRange12 const idx_range_eval_periodic_rtheta2(idx_range_eval_r2, idx_range_theta2);
    FieldMem<Patch2::Coord12, Patch2::IdxRange12> eval_points_2_alloc(
            idx_range_eval_periodic_rtheta2);
    Field<Patch2::Coord12, Patch2::IdxRange12> eval_points_2 = get_field(eval_points_2_alloc);

    set_eval_points_2D<Patch1::Grid1, Patch1::Grid2>(eval_points_1);
    set_eval_points_2D<Patch2::Grid1, Patch2::Grid2>(eval_points_2);

    // --- collection
    MultipatchField<CoordConstFieldOnPatch, Patch1, Patch2> const
            eval_points(get_const_field(eval_points_1), get_const_field(eval_points_2));


    // Evaluated functions
    // --- patch 1
    DFieldMem<Patch1::IdxRange12> eval_function_1_alloc(idx_range_eval_periodic_rtheta1);
    DField<Patch1::IdxRange12> eval_function_1 = get_field(eval_function_1_alloc);

    // --- patch 2
    DFieldMem<Patch2::IdxRange12> eval_function_2_alloc(idx_range_eval_periodic_rtheta2);
    DField<Patch2::IdxRange12> eval_function_2 = get_field(eval_function_2_alloc);

    // --- collection
    MultipatchField<DFieldOnPatch, Patch1, Patch2> const
            eval_functions(eval_function_1, eval_function_2);


    // The builder collection owned by the interpolation must give the spline coefficients
    // built in the fixture by hand-built builders (f = r sin(theta) and f = r^2 sin(theta)).
    MultipatchFieldMem<DFieldMemOnPatch, Patch1, Patch2> values_alloc(all_idx_ranges);
    initialise_function<Patch1>(values_alloc.get<Patch1>(), 1);
    initialise_function<Patch2>(values_alloc.get<Patch2>(), 2);
    MultipatchField<DConstFieldOnPatch, Patch1, Patch2> const values
            = get_const_field(values_alloc);

    MultipatchFieldMem<SplineCoeffMemOnPatch_2D, Patch1, Patch2> coeffs_alloc(
            interpolation.get_builder().batched_basis_idx_ranges(all_idx_ranges));
    MultipatchField<SplineCoeffOnPatch_2D, Patch1, Patch2> const coeffs = get_field(coeffs_alloc);
    interpolation.get_builder()(coeffs, values);

    EXPECT_LE(
            max_abs_difference(get_const_field(coeffs_alloc.get<Patch1>()), splines.get<Patch1>()),
            1e-14);
    EXPECT_LE(
            max_abs_difference(get_const_field(coeffs_alloc.get<Patch2>()), splines.get<Patch2>()),
            1e-14);

    // Evaluate the functions at the evaluation points.
    evaluators(eval_functions, eval_points, splines);

    // Expected functions
    // --- patch 1
    DFieldMem<Patch1::IdxRange12> expected_function_1_alloc(idx_range_eval_periodic_rtheta1);
    DField<Patch1::IdxRange12> expected_function_1 = get_field(expected_function_1_alloc);

    // --- patch 2
    DFieldMem<Patch2::IdxRange12> expected_function_2_alloc(idx_range_eval_periodic_rtheta2);
    DField<Patch2::IdxRange12> expected_function_2 = get_field(expected_function_2_alloc);

    evaluator_1(
            expected_function_1,
            get_const_field(eval_points_1),
            get_const_field(function_1_coef));
    evaluator_2(
            expected_function_2,
            get_const_field(eval_points_2),
            get_const_field(function_2_coef));


    // Compare the evaluated functions with the expected functions.
    // --- put the functions on host
    auto const expected_function_1_host = ddc::create_mirror_and_copy(expected_function_1);
    auto const expected_function_2_host = ddc::create_mirror_and_copy(expected_function_2);
    auto const eval_function_1_host = ddc::create_mirror_and_copy(eval_function_1);
    auto const eval_function_2_host = ddc::create_mirror_and_copy(eval_function_2);

    // --- check errors
    ddc::host_for_each(idx_range_eval_periodic_rtheta1, [&](typename Patch1::Idx12 const idx) {
        EXPECT_NEAR(eval_function_1_host(idx), expected_function_1_host(idx), 1e-12);
    });
    ddc::host_for_each(idx_range_eval_periodic_rtheta2, [&](typename Patch2::Idx12 const idx) {
        EXPECT_NEAR(eval_function_2_host(idx), expected_function_2_host(idx), 1e-12);
    });
}


/* -----------------------------------------------------------------------------------------------
    Test deriv() for a single coordinate on host.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionTest, HostDerivativesOnSingleCoord)
{
    // Compare the evaluated functions with the expected functions.
    // Evaluation points
    // --- patch 1
    typename Patch1::Coord12 const eval_coord_P1_on_patch_1(0.5, 0.0);
    typename Patch1::Coord12 const eval_coord_P1_on_patch_2(1.1, M_PI);

    // --- patch 2
    typename Patch2::Coord12 const eval_coord_P2_on_patch_1(0.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 const eval_coord_P2_on_patch_2(1.9, 0.0);


    auto const function_1_coef_host = ddc::create_mirror_and_copy(function_1_coef);
    auto const function_2_coef_host = ddc::create_mirror_and_copy(function_2_coef);
    MultipatchField<ConstSplineCoeffOnPatch_2D_host, Patch1, Patch2> const
            splines_host(get_field(function_1_coef_host), get_field(function_2_coef_host));

    ConstSplineCoeffOnPatch_2D_host<Patch1> const spline_1_host(
            get_const_field(function_1_coef_host));
    ConstSplineCoeffOnPatch_2D_host<Patch2> const spline_2_host(
            get_const_field(function_2_coef_host));

    Kokkos::View<Coord<R, Theta>*, Kokkos::HostSpace> coords_host("coords_host", 4);

    coords_host(0) = eval_coord_P1_on_patch_1;
    coords_host(1) = eval_coord_P1_on_patch_2;
    coords_host(2) = eval_coord_P2_on_patch_1;
    coords_host(3) = eval_coord_P2_on_patch_2;

    // Check deriv() on host.
    test_deriv<ConstSplineCoeffOnPatch_2D_host>(
            evaluators_host,
            evaluator_1_host,
            evaluator_2_host,
            coords_host,
            splines_host,
            spline_1_host,
            spline_2_host);
}


/* -----------------------------------------------------------------------------------------------
    Test deriv for a single coordinate on device.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionTest, DeviceDerivativesOnSingleCoord)
{
    // Compare the evaluated functions with the expected functions.
    // Evaluation points
    // --- patch 1
    typename Patch1::Coord12 const eval_coord_P1_on_patch_1(0.5, 0.0);
    typename Patch1::Coord12 const eval_coord_P1_on_patch_2(1.1, M_PI);

    // --- patch 2
    typename Patch2::Coord12 const eval_coord_P2_on_patch_1(0.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 const eval_coord_P2_on_patch_2(1.9, 0.0);

    Kokkos::View<Coord<R, Theta>*> coords("coords", 4);
    Kokkos::View<Coord<R, Theta>*, Kokkos::HostSpace> coords_host("coords_host", 4);

    coords_host(0) = eval_coord_P1_on_patch_1;
    coords_host(1) = eval_coord_P1_on_patch_2;
    coords_host(2) = eval_coord_P2_on_patch_1;
    coords_host(3) = eval_coord_P2_on_patch_2;

    Kokkos::deep_copy(coords, coords_host);

    // Check deriv() on device.
    test_deriv<ConstSplineCoeffOnPatch_2D>(
            evaluators,
            evaluator_1,
            evaluator_2,
            coords,
            splines,
            get_const_field(function_1_coef),
            get_const_field(function_2_coef));
}


/* -----------------------------------------------------------------------------------------------
    Test deriv() for single coordinates outside of the domain. The operators are supposed to
    return an error.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionTest, DerivativesOnSingleCoordDeathTest)
{
    // Evaluation points
    typename Patch1::Coord12 const eval_coord_P1_outside(2.5, 0.0);
    typename Patch2::Coord12 const eval_coord_P2_outside(2.5, 0.0);

    Idx<ddc::Deriv<R>> d_r(1);
    Idx<ddc::Deriv<Theta>> d_theta(1);
    Idx<ddc::Deriv<R>, ddc::Deriv<Theta>> d_rtheta(1, 1);

    // Test the Kokkos::abort only on host.
    if constexpr (std::is_same_v<
                          Kokkos::DefaultExecutionSpace::memory_space,
                          Kokkos::DefaultHostExecutionSpace::memory_space>) {
        // --- derivativative 1
        EXPECT_DEATH(
                evaluators.deriv(d_r, eval_coord_P1_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");

        EXPECT_DEATH(
                evaluators.deriv(d_r, eval_coord_P2_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");

        // --- derivativative 2
        EXPECT_DEATH(
                evaluators.deriv(d_theta, eval_coord_P1_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");

        EXPECT_DEATH(
                evaluators.deriv(d_theta, eval_coord_P2_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");

        // --- cross-derivativative
        EXPECT_DEATH(
                evaluators.deriv(d_rtheta, eval_coord_P1_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");

        EXPECT_DEATH(
                evaluators.deriv(d_rtheta, eval_coord_P2_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");
    }
}


/* -----------------------------------------------------------------------------------------------
    Test deriv() for fields of coordinates called from host and stored on device.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionTest, DerivativesOnCoordField)
{
    // Evaluation points
    // --- patch 1
    FieldMem<Patch1::Coord12, Patch1::IdxRange12> eval_points_1_alloc(idx_range_eval_rtheta1);
    Field<Patch1::Coord12, Patch1::IdxRange12> eval_points_1 = get_field(eval_points_1_alloc);

    // --- patch 1
    FieldMem<Patch2::Coord12, Patch2::IdxRange12> eval_points_2_alloc(idx_range_eval_rtheta2);
    Field<Patch2::Coord12, Patch2::IdxRange12> eval_points_2 = get_field(eval_points_2_alloc);

    set_eval_points_2D<Patch1::Grid1, Patch1::Grid2>(eval_points_1);
    set_eval_points_2D<Patch2::Grid1, Patch2::Grid2>(eval_points_2);

    ConstField<Patch1::Coord12, Patch1::IdxRange12> const_eval_points_1
            = get_const_field(eval_points_1);
    ConstField<Patch2::Coord12, Patch2::IdxRange12> const_eval_points_2
            = get_const_field(eval_points_2);

    // --- collection
    MultipatchField<CoordConstFieldOnPatch, Patch1, Patch2> const
            eval_points(const_eval_points_1, const_eval_points_2);


    // Evaluated functions
    // --- derivativative 1
    DFieldMem<Patch1::IdxRange12> eval_derivs_1_patch_1_alloc(idx_range_eval_rtheta1);
    DField<Patch1::IdxRange12> eval_derivs_1_patch_1 = get_field(eval_derivs_1_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> eval_derivs_1_patch_2_alloc(idx_range_eval_rtheta2);
    DField<Patch2::IdxRange12> eval_derivs_1_patch_2 = get_field(eval_derivs_1_patch_2_alloc);

    // --- derivativative 2
    DFieldMem<Patch1::IdxRange12> eval_derivs_2_patch_1_alloc(idx_range_eval_rtheta1);
    DField<Patch1::IdxRange12> eval_derivs_2_patch_1 = get_field(eval_derivs_2_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> eval_derivs_2_patch_2_alloc(idx_range_eval_rtheta2);
    DField<Patch2::IdxRange12> eval_derivs_2_patch_2 = get_field(eval_derivs_2_patch_2_alloc);

    // --- cross-derivatives
    DFieldMem<Patch1::IdxRange12> eval_derivs_12_patch_1_alloc(idx_range_eval_rtheta1);
    DField<Patch1::IdxRange12> eval_derivs_12_patch_1 = get_field(eval_derivs_12_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> eval_derivs_12_patch_2_alloc(idx_range_eval_rtheta2);
    DField<Patch2::IdxRange12> eval_derivs_12_patch_2 = get_field(eval_derivs_12_patch_2_alloc);

    // ----- collections
    MultipatchField<DFieldOnPatch, Patch1, Patch2> const
            eval_derivs_1(eval_derivs_1_patch_1, eval_derivs_1_patch_2);
    MultipatchField<DFieldOnPatch, Patch1, Patch2> const
            eval_derivs_2(eval_derivs_2_patch_1, eval_derivs_2_patch_2);
    MultipatchField<DFieldOnPatch, Patch1, Patch2> const
            eval_derivs_12(eval_derivs_12_patch_1, eval_derivs_12_patch_2);


    // Derivative orders (arbitrary between 1 and 2)
    Idx<ddc::Deriv<R>> d_r(1);
    Idx<ddc::Deriv<Theta>> d_theta(2);
    Idx<ddc::Deriv<R>, ddc::Deriv<Theta>> d_r_d_theta(2, 1);

    // Evaluate the functions at the evaluation points.
    evaluators.deriv(d_r, eval_derivs_1, eval_points, splines);
    evaluators.deriv(d_theta, eval_derivs_2, eval_points, splines);
    evaluators.deriv(d_r_d_theta, eval_derivs_12, eval_points, splines);


    // Expected functions
    // --- derivative 1
    DFieldMem<Patch1::IdxRange12> expected_derivs_1_patch_1_alloc(idx_range_eval_rtheta1);
    DField<Patch1::IdxRange12> expected_derivs_1_patch_1
            = get_field(expected_derivs_1_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> expected_derivs_1_patch_2_alloc(idx_range_eval_rtheta2);
    DField<Patch2::IdxRange12> expected_derivs_1_patch_2
            = get_field(expected_derivs_1_patch_2_alloc);

    // --- derivative 2
    DFieldMem<Patch1::IdxRange12> expected_derivs_2_patch_1_alloc(idx_range_eval_rtheta1);
    DField<Patch1::IdxRange12> expected_derivs_2_patch_1
            = get_field(expected_derivs_2_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> expected_derivs_2_patch_2_alloc(idx_range_eval_rtheta2);
    DField<Patch2::IdxRange12> expected_derivs_2_patch_2
            = get_field(expected_derivs_2_patch_2_alloc);

    // --- cross-derivative
    DFieldMem<Patch1::IdxRange12> expected_derivs_12_patch_1_alloc(idx_range_eval_rtheta1);
    DField<Patch1::IdxRange12> expected_derivs_12_patch_1
            = get_field(expected_derivs_12_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> expected_derivs_12_patch_2_alloc(idx_range_eval_rtheta2);
    DField<Patch2::IdxRange12> expected_derivs_12_patch_2
            = get_field(expected_derivs_12_patch_2_alloc);


    DConstField<typename Patch1::IdxRangeBS12> const const_function_1_coef
            = get_const_field(function_1_coef);
    DConstField<typename Patch2::IdxRangeBS12> const const_function_2_coef
            = get_const_field(function_2_coef);


    evaluator_1.deriv(d_r, expected_derivs_1_patch_1, const_eval_points_1, const_function_1_coef);
    evaluator_2.deriv(d_r, expected_derivs_1_patch_2, const_eval_points_2, const_function_2_coef);

    evaluator_1
            .deriv(d_theta, expected_derivs_2_patch_1, const_eval_points_1, const_function_1_coef);
    evaluator_2
            .deriv(d_theta, expected_derivs_2_patch_2, const_eval_points_2, const_function_2_coef);

    evaluator_1
            .deriv(d_r_d_theta,
                   expected_derivs_12_patch_1,
                   const_eval_points_1,
                   const_function_1_coef);
    evaluator_2
            .deriv(d_r_d_theta,
                   expected_derivs_12_patch_2,
                   const_eval_points_2,
                   const_function_2_coef);


    // Compare the evaluated derivatives with the expected derivatives.
    // --- put the functions on host
    auto const expected_derivs_1_patch_1_host
            = ddc::create_mirror_and_copy(expected_derivs_1_patch_1);
    auto const expected_derivs_1_patch_2_host
            = ddc::create_mirror_and_copy(expected_derivs_1_patch_2);
    auto const eval_derivs_1_patch_1_host = ddc::create_mirror_and_copy(eval_derivs_1_patch_1);
    auto const eval_derivs_1_patch_2_host = ddc::create_mirror_and_copy(eval_derivs_1_patch_2);

    auto const expected_derivs_2_patch_1_host
            = ddc::create_mirror_and_copy(expected_derivs_2_patch_1);
    auto const expected_derivs_2_patch_2_host
            = ddc::create_mirror_and_copy(expected_derivs_2_patch_2);
    auto const eval_derivs_2_patch_1_host = ddc::create_mirror_and_copy(eval_derivs_2_patch_1);
    auto const eval_derivs_2_patch_2_host = ddc::create_mirror_and_copy(eval_derivs_2_patch_2);

    auto const expected_derivs_12_patch_1_host
            = ddc::create_mirror_and_copy(expected_derivs_12_patch_1);
    auto const expected_derivs_12_patch_2_host
            = ddc::create_mirror_and_copy(expected_derivs_12_patch_2);
    auto const eval_derivs_12_patch_1_host = ddc::create_mirror_and_copy(eval_derivs_12_patch_1);
    auto const eval_derivs_12_patch_2_host = ddc::create_mirror_and_copy(eval_derivs_12_patch_2);

    // --- check errors
    ddc::host_for_each(idx_range_eval_rtheta1, [&](typename Patch1::Idx12 const idx) {
        EXPECT_NEAR(expected_derivs_1_patch_1_host(idx), eval_derivs_1_patch_1_host(idx), 1e-11);
        EXPECT_NEAR(expected_derivs_2_patch_1_host(idx), eval_derivs_2_patch_1_host(idx), 1e-11);
        EXPECT_NEAR(expected_derivs_12_patch_1_host(idx), eval_derivs_12_patch_1_host(idx), 1e-10);
    });
    ddc::host_for_each(idx_range_eval_rtheta2, [&](typename Patch2::Idx12 const idx) {
        EXPECT_NEAR(expected_derivs_1_patch_2_host(idx), eval_derivs_1_patch_2_host(idx), 1e-11);
        EXPECT_NEAR(expected_derivs_2_patch_2_host(idx), eval_derivs_2_patch_2_host(idx), 1e-11);
        EXPECT_NEAR(expected_derivs_12_patch_2_host(idx), eval_derivs_12_patch_2_host(idx), 1e-10);
    });
}


/* -----------------------------------------------------------------------------------------------
    Test integrate() on host.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionTest, HostIntegrateOnCoordField)
{
    auto function_1_coef_host = ddc::create_mirror_and_copy(function_1_coef);
    auto function_2_coef_host = ddc::create_mirror_and_copy(function_2_coef);
    MultipatchField<ConstSplineCoeffOnPatch_2D_host, Patch1, Patch2> const splines_host(
            get_const_field(function_1_coef_host),
            get_const_field(function_2_coef_host));

    // Evaluate the functions at the evaluation points.
    Kokkos::View<double[2], HostExecSpace> eval_integrals("eval_integrals", 2);
    evaluators_host.integrate(eval_integrals, splines_host);

    // Expected functions
    double const expected_integral_1 = 0.0;
    double const expected_integral_2 = 0.0;

    // Compare the evaluated functions with the expected functions.
    EXPECT_NEAR(eval_integrals(0), expected_integral_1, 1e-12);
    EXPECT_NEAR(eval_integrals(1), expected_integral_2, 1e-12);

    /*
        REMARK: Cannot compare to SplineEvaluator2D as long as the .integrate()
            function is only defined on batch domain.
    */
}


/* -----------------------------------------------------------------------------------------------
    Test operator() for DerivField and fields of coordinates called from host but stored on device.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionTest, EvaluateDerivFieldOnCoordField)
{
    // Evaluation points
    // --- patch 1
    Patch1::IdxRange12 const idx_range_eval_periodic_rtheta1(idx_range_eval_r1, idx_range_theta1);
    FieldMem<Patch1::Coord12, Patch1::IdxRange12> eval_points_1_alloc(
            idx_range_eval_periodic_rtheta1);
    Field<Patch1::Coord12, Patch1::IdxRange12> eval_points_1 = get_field(eval_points_1_alloc);

    // --- patch 1
    Patch2::IdxRange12 const idx_range_eval_periodic_rtheta2(idx_range_eval_r2, idx_range_theta2);
    FieldMem<Patch2::Coord12, Patch2::IdxRange12> eval_points_2_alloc(
            idx_range_eval_periodic_rtheta2);
    Field<Patch2::Coord12, Patch2::IdxRange12> eval_points_2 = get_field(eval_points_2_alloc);

    set_eval_points_2D<Patch1::Grid1, Patch1::Grid2>(eval_points_1);
    set_eval_points_2D<Patch2::Grid1, Patch2::Grid2>(eval_points_2);

    // --- collection
    MultipatchField<CoordConstFieldOnPatch, Patch1, Patch2> const
            eval_points(get_const_field(eval_points_1), get_const_field(eval_points_2));


    // Evaluated functions
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2>
            idx_range_evals(idx_range_eval_periodic_rtheta1, idx_range_eval_periodic_rtheta2);

    MultipatchFieldMem<DerivFieldMemOnPatch, Patch1, Patch2> eval_functions_and_derivs_alloc(
            idx_range_evals);
    MultipatchField<DerivFieldOnPatch, Patch1, Patch2> eval_functions_and_derivs(
            eval_functions_and_derivs_alloc);

    // Evaluate the functions at the evaluation points.
    evaluators(eval_functions_and_derivs, eval_points, splines);

    // Expected functions
    // --- patch 1
    DFieldMem<Patch1::IdxRange12> expected_function_1_alloc(idx_range_eval_periodic_rtheta1);
    DField<Patch1::IdxRange12> expected_function_1 = get_field(expected_function_1_alloc);

    // --- patch 2
    DFieldMem<Patch2::IdxRange12> expected_function_2_alloc(idx_range_eval_periodic_rtheta2);
    DField<Patch2::IdxRange12> expected_function_2 = get_field(expected_function_2_alloc);

    evaluator_1(
            expected_function_1,
            get_const_field(eval_points_1),
            get_const_field(function_1_coef));
    evaluator_2(
            expected_function_2,
            get_const_field(eval_points_2),
            get_const_field(function_2_coef));


    // Compare the evaluated functions with the expected functions.
    // --- put the functions on host
    auto const expected_function_1_host = ddc::create_mirror_and_copy(expected_function_1);
    auto const expected_function_2_host = ddc::create_mirror_and_copy(expected_function_2);

    // --- define DField for ddc::create_mirror_and_copy().
    // --- patch 1
    DFieldMem<Patch1::IdxRange12> eval_function_1_alloc(idx_range_eval_periodic_rtheta1);
    DField<Patch1::IdxRange12> eval_function_1 = get_field(eval_function_1_alloc);

    // --- patch 2
    DFieldMem<Patch2::IdxRange12> eval_function_2_alloc(idx_range_eval_periodic_rtheta2);
    DField<Patch2::IdxRange12> eval_function_2 = get_field(eval_function_2_alloc);

    ddc::parallel_deepcopy(
            eval_function_1,
            eval_functions_and_derivs.template get<Patch1>().get_values_field());
    ddc::parallel_deepcopy(
            eval_function_2,
            eval_functions_and_derivs.template get<Patch2>().get_values_field());

    auto const eval_function_1_host = ddc::create_mirror_and_copy(eval_function_1);
    auto const eval_function_2_host = ddc::create_mirror_and_copy(eval_function_2);

    // --- check errors
    ddc::host_for_each(idx_range_eval_periodic_rtheta1, [&](typename Patch1::Idx12 const idx) {
        EXPECT_NEAR(eval_function_1_host(idx), expected_function_1_host(idx), 1e-12);
    });
    ddc::host_for_each(idx_range_eval_periodic_rtheta2, [&](typename Patch2::Idx12 const idx) {
        EXPECT_NEAR(eval_function_2_host(idx), expected_function_2_host(idx), 1e-12);
    });
}


/* -----------------------------------------------------------------------------------------------
    Test deriv() for DerivField and fields of coordinates called from host and stored on device.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionTest, DerivativesDerivFieldOnCoordField)
{
    // Evaluation points
    // --- patch 1
    FieldMem<Patch1::Coord12, Patch1::IdxRange12> eval_points_1_alloc(idx_range_eval_rtheta1);
    Field<Patch1::Coord12, Patch1::IdxRange12> eval_points_1 = get_field(eval_points_1_alloc);

    // --- patch 1
    FieldMem<Patch2::Coord12, Patch2::IdxRange12> eval_points_2_alloc(idx_range_eval_rtheta2);
    Field<Patch2::Coord12, Patch2::IdxRange12> eval_points_2 = get_field(eval_points_2_alloc);

    set_eval_points_2D<Patch1::Grid1, Patch1::Grid2>(eval_points_1);
    set_eval_points_2D<Patch2::Grid1, Patch2::Grid2>(eval_points_2);

    // --- collection
    MultipatchField<CoordConstFieldOnPatch, Patch1, Patch2> const
            eval_points(get_const_field(eval_points_1), get_const_field(eval_points_2));


    // Evaluated functions
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2>
            idx_range_evals(idx_range_eval_rtheta1, idx_range_eval_rtheta2);

    MultipatchFieldMem<DerivFieldMemOnPatch, Patch1, Patch2> eval_functions_and_derivs_alloc(
            idx_range_evals);
    MultipatchField<DerivFieldOnPatch, Patch1, Patch2> eval_functions_and_derivs(
            eval_functions_and_derivs_alloc);

    // Derivative orders.
    Idx<ddc::Deriv<R>> d_r(1);
    Idx<ddc::Deriv<Theta>> d_theta(0);
    Idx<ddc::Deriv<R>, ddc::Deriv<Theta>> d_r_d_theta(1, 1);

    // Evaluate the functions at the evaluation points.
    evaluators.deriv(d_r, eval_functions_and_derivs, eval_points, splines);
    evaluators.deriv(d_theta, eval_functions_and_derivs, eval_points, splines);
    evaluators.deriv(d_r_d_theta, eval_functions_and_derivs, eval_points, splines);

    DConstField<typename Patch1::IdxRangeBS12> const const_spline_1
            = get_const_field(function_1_coef);
    DConstField<typename Patch2::IdxRangeBS12> const const_spline_2
            = get_const_field(function_2_coef);


    // Check derivatives.
    // --- derivatives 1
    check_deriv_match<Patch1, typename Patch1::Grid1>(
            d_r,
            idx_range_eval_rtheta1,
            eval_functions_and_derivs.template get<Patch1>(),
            evaluator_1,
            const_spline_1);
    check_deriv_match<Patch2, typename Patch2::Grid1>(
            d_r,
            idx_range_eval_rtheta2,
            eval_functions_and_derivs.template get<Patch2>(),
            evaluator_2,
            const_spline_2);

    // --- derivatives 2
    check_deriv_match<Patch1, typename Patch1::Grid2>(
            d_theta,
            idx_range_eval_rtheta1,
            eval_functions_and_derivs.template get<Patch1>(),
            evaluator_1,
            const_spline_1);
    check_deriv_match<Patch2, typename Patch2::Grid2>(
            d_theta,
            idx_range_eval_rtheta2,
            eval_functions_and_derivs.template get<Patch2>(),
            evaluator_2,
            const_spline_2);

    // --- cross-derivatives
    check_cross_deriv_match<Patch1>(
            d_r_d_theta,
            idx_range_eval_rtheta1,
            eval_functions_and_derivs.template get<Patch1>(),
            evaluator_1,
            const_spline_1);
    check_cross_deriv_match<Patch2>(
            d_r_d_theta,
            idx_range_eval_rtheta2,
            eval_functions_and_derivs.template get<Patch2>(),
            evaluator_2,
            const_spline_2);
}


/* -----------------------------------------------------------------------------------------------
    Test operator() for single coordinates with non-spline (Lagrange) evaluators.
    The coordinates are located on different patches, outside the periodic domain and
    outside the global domain.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionLagrangeTest, EvaluateOnSingleCoord)
{
    std::array<Coord<R, Theta>, 5> const eval_coords_array {
            Coord<R, Theta>(0.5, 0.0),
            Coord<R, Theta>(1.1, M_PI),
            Coord<R, Theta>(0.5, 3 / 2. * M_PI),
            Coord<R, Theta>(1.9, 4 * M_PI),
            Coord<R, Theta>(2.3, -M_PI)};

    Kokkos::View<Coord<R, Theta>*, Kokkos::DefaultExecutionSpace::memory_space>
            coords("coords", eval_coords_array.size());
    Kokkos::View<Coord<R, Theta>*, Kokkos::HostSpace>
            coords_host("coords_host", eval_coords_array.size());
    for (std::size_t i(0); i < eval_coords_array.size(); ++i) {
        coords_host(i) = eval_coords_array[i];
    }
    Kokkos::deep_copy(coords, coords_host);

    check_lagrange_on_single_coords(evaluators, evaluator_1, evaluator_2, coords, coeffs);
}


/* -----------------------------------------------------------------------------------------------
    Test operator() and deriv() for fields of coordinates with non-spline (Lagrange) evaluators.
    --------------------------------------------------------------------------------------------*/
TEST_F(InterpolationEvaluator2DCollectionLagrangeTest, EvaluateOnCoordField)
{
    // Evaluation points
    Patch1::IdxRange12 const idx_range_eval_1(
            Patch1::IdxRange1(idx_range_rtheta1).remove_last(Patch1::IdxStep1(1)),
            Patch1::IdxRange2(idx_range_rtheta1));
    Patch2::IdxRange12 const idx_range_eval_2(
            Patch2::IdxRange1(idx_range_rtheta2).remove_last(Patch2::IdxStep1(1)),
            Patch2::IdxRange2(idx_range_rtheta2));

    MultipatchType<IdxRangeOnPatch, Patch1, Patch2> idx_range(idx_range_eval_1, idx_range_eval_2);

    MultipatchFieldMem<CoordFieldMemOnPatch, Patch1, Patch2> eval_points_alloc(idx_range);
    set_eval_points_2D<Patch1::Grid1, Patch1::Grid2>(eval_points_alloc.get<Patch1>());
    set_eval_points_2D<Patch2::Grid1, Patch2::Grid2>(eval_points_alloc.get<Patch2>());

    MultipatchField<CoordConstFieldOnPatch, Patch1, Patch2> const eval_points
            = get_const_field(eval_points_alloc);

    // Values
    check_lagrange_on_coord_field(
            Idx<>(),
            evaluators,
            evaluator_1,
            evaluator_2,
            eval_points,
            coeffs);

    // Derivatives
    check_lagrange_on_coord_field(
            Idx<ddc::Deriv<R>>(1),
            evaluators,
            evaluator_1,
            evaluator_2,
            eval_points,
            coeffs);
    check_lagrange_on_coord_field(
            Idx<ddc::Deriv<Theta>>(2),
            evaluators,
            evaluator_1,
            evaluator_2,
            eval_points,
            coeffs);
    check_lagrange_on_coord_field(
            Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1),
            evaluators,
            evaluator_1,
            evaluator_2,
            eval_points,
            coeffs);
}
