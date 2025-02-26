// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_non_uniform.hpp"
#include "constant_extrapolation_rules_onion.hpp"
#include "mesh_builder.hpp"
#include "multipatch_field.hpp"
#include "multipatch_spline_evaluator_2d.hpp"
#include "multipatch_type.hpp"
#include "onion_patch_locator.hpp"
#include "physical_geometry.hpp"
#include "spline_testing_tools.hpp"
#include "types.hpp"


namespace {
using DeviceMultipatchSplineRThetaEvaluator = MultipatchSplineEvaluator2D<
        DeviceExecSpace,
        typename DeviceExecSpace::memory_space,
        BSplines1OnPatch,
        BSplines2OnPatch,
        Grid1OnPatch,
        Grid2OnPatch,
        ConstantExtrapolationRuleOnion<PatchLocator<DeviceExecSpace>>,
        DFieldOnPatch,
        PatchLocator<DeviceExecSpace>,
        Patch1,
        Patch2>;

using HostMultipatchSplineRThetaEvaluator = MultipatchSplineEvaluator2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplines1OnPatch,
        BSplines2OnPatch,
        Grid1OnPatch,
        Grid2OnPatch,
        ConstantExtrapolationRuleOnion<PatchLocator<HostExecSpace>>,
        DFieldOnPatch_host,
        PatchLocator<HostExecSpace>,
        Patch1,
        Patch2>;


class MultipatchSplineEvaluatorTest : public MultipatchSplineOnionShapeTest
{
protected:
    // Index ranges
    Patch1::IdxRange1 const idx_range_r1;
    Patch1::IdxRange2 const idx_range_theta1;

    Patch2::IdxRange1 const idx_range_r2;
    Patch2::IdxRange2 const idx_range_theta2;

    // Single extrapolation rules
    ddc::ConstantExtrapolationRule<R, Theta> const bc_r_min_1;
    ddc::ConstantExtrapolationRule<R, Theta> const bc_r_max_1;
    ddc::ConstantExtrapolationRule<R, Theta> const bc_r_min_2;
    ddc::ConstantExtrapolationRule<R, Theta> const bc_r_max_2;
    ddc::PeriodicExtrapolationRule<Theta> const bc_theta;

    // Evaluators on single patch
    SplineRThetaEvaluator<1, DeviceExecSpace> const evaluator_1;
    SplineRThetaEvaluator<2, DeviceExecSpace> const evaluator_2;

    SplineRThetaEvaluator<1, HostExecSpace> const evaluator_1_host;
    SplineRThetaEvaluator<2, HostExecSpace> const evaluator_2_host;

    // Spline representations
    DField<IdxRange<BSplinesR<1>, BSplinesTheta<1>>> const function_1_coef;
    DField<IdxRange<BSplinesR<2>, BSplinesTheta<2>>> const function_2_coef;

    LogicalToPhysicalMapping const to_physical_mapping;
    PhysicalToLogicalMapping const to_logical_mapping;

public:
    MultipatchSplineEvaluatorTest()
        // Idx ranges
        : idx_range_r1(SplineInterpPointsR<1>::get_domain<GridR<1>>())
        , idx_range_theta1(SplineInterpPointsTheta<1>::get_domain<GridTheta<1>>())
        , idx_range_r2(SplineInterpPointsR<2>::get_domain<GridR<2>>())
        , idx_range_theta2(SplineInterpPointsTheta<2>::get_domain<GridTheta<2>>())
        , bc_r_min_1(ddc::coordinate(idx_range_r1.front()))
        , bc_r_max_1(ddc::coordinate(idx_range_r1.back()))
        , bc_r_min_2(ddc::coordinate(idx_range_r2.front()))
        , bc_r_max_2(ddc::coordinate(idx_range_r2.back()))
        , bc_theta()
        // Local spline evaluators on device and host
        , evaluator_1(bc_r_min_1, bc_r_max_1, bc_theta, bc_theta)
        , evaluator_2(bc_r_min_2, bc_r_max_2, bc_theta, bc_theta)
        , evaluator_1_host(bc_r_min_1, bc_r_max_1, bc_theta, bc_theta)
        , evaluator_2_host(bc_r_min_2, bc_r_max_2, bc_theta, bc_theta)
        // Local splines on device and host
        , function_1_coef(get_field(function_1_coef_alloc))
        , function_2_coef(get_field(function_2_coef_alloc))
        , to_physical_mapping()
        , to_logical_mapping()
    {
    }


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

        IdxRange<GridR, GridTheta> const reduced_idx_range = get_idx_range(eval_points);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                reduced_idx_range,
                KOKKOS_LAMBDA(Idx<GridR, GridTheta> const idx) {
                    CoordRTheta const coord = ddc::coordinate(idx);
                    CoordR const delta_r = distance_at_right(Idx<GridR>(idx));
                    CoordTheta const delta_theta = distance_at_right(Idx<GridTheta>(idx));
                    eval_points(idx) = coord + CoordRTheta(delta_r / 2., delta_theta / 2.);
                });
    }
};


// Checking functions ----------------------------------------------------------------------------
// One function per method of MultipatchSplineEvaluator tested.
// Test operator() ...............................................................................
void test_operator_assignement(
        DeviceMultipatchSplineRThetaEvaluator const& evaluators,
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
                err = Kokkos::max(abs(eval_function - expected_function), err);
                err = 0;
            },
            Kokkos::Max<double>(max_error));
    EXPECT_LE(max_error, 1e-15);
}


// Test deriv_dim_1() ............................................................................
void test_deriv_dim_1(
        DeviceMultipatchSplineRThetaEvaluator const& evaluators,
        SplineRThetaEvaluator<1, DeviceExecSpace> const& single_evaluator_1,
        SplineRThetaEvaluator<2, DeviceExecSpace> const& single_evaluator_2,
        Kokkos::View<Coord<R, Theta>*> const& eval_coords,
        MultipatchField<ConstSplineCoeffOnPatch_2D, Patch1, Patch2> const& splines,
        ConstSplineCoeffOnPatch_2D<Patch1> const& spline_patch_1,
        ConstSplineCoeffOnPatch_2D<Patch2> const& spline_patch_2)
{
    double max_error = 0;
    Kokkos::parallel_reduce(
            eval_coords.extent(0),
            KOKKOS_LAMBDA(int i, double& err) {
                double const eval_function = evaluators.deriv_dim_1(eval_coords(i), splines);
                double expected_function;
                if (Coord<R>(eval_coords(i)) < ddc::discrete_space<BSplinesR<1>>().rmax()) {
                    expected_function
                            = single_evaluator_1.deriv_dim_1(eval_coords(i), spline_patch_1);
                } else {
                    expected_function
                            = single_evaluator_2.deriv_dim_1(eval_coords(i), spline_patch_2);
                }
                err = Kokkos::max(abs(eval_function - expected_function), err);
            },
            Kokkos::Max<double>(max_error));
    EXPECT_LE(max_error, 1e-15);
};

// Test deriv_dim_2() ............................................................................
void test_deriv_dim_2(
        DeviceMultipatchSplineRThetaEvaluator const& evaluators,
        SplineRThetaEvaluator<1, DeviceExecSpace> const& single_evaluator_1,
        SplineRThetaEvaluator<2, DeviceExecSpace> const& single_evaluator_2,
        Kokkos::View<Coord<R, Theta>*> const& eval_coords,
        MultipatchField<ConstSplineCoeffOnPatch_2D, Patch1, Patch2> const& splines,
        ConstSplineCoeffOnPatch_2D<Patch1> const& spline_patch_1,
        ConstSplineCoeffOnPatch_2D<Patch2> const& spline_patch_2)
{
    double max_error = 0;
    Kokkos::parallel_reduce(
            eval_coords.extent(0),
            KOKKOS_LAMBDA(int i, double& err) {
                double const eval_function = evaluators.deriv_dim_2(eval_coords(i), splines);
                double expected_function;
                if (Coord<R>(eval_coords(i)) < ddc::discrete_space<BSplinesR<1>>().rmax()) {
                    expected_function
                            = single_evaluator_1.deriv_dim_2(eval_coords(i), spline_patch_1);
                } else {
                    expected_function
                            = single_evaluator_2.deriv_dim_2(eval_coords(i), spline_patch_2);
                }
                err = Kokkos::max(abs(eval_function - expected_function), err);
            },
            Kokkos::Max<double>(max_error));
    EXPECT_LE(max_error, 1e-15);
};

// Test deriv_1_and_2() ..........................................................................
void test_deriv_1_and_2(
        DeviceMultipatchSplineRThetaEvaluator const& evaluators,
        SplineRThetaEvaluator<1, DeviceExecSpace> const& single_evaluator_1,
        SplineRThetaEvaluator<2, DeviceExecSpace> const& single_evaluator_2,
        Kokkos::View<Coord<R, Theta>*> const& eval_coords,
        MultipatchField<ConstSplineCoeffOnPatch_2D, Patch1, Patch2> const& splines,
        ConstSplineCoeffOnPatch_2D<Patch1> const& spline_patch_1,
        ConstSplineCoeffOnPatch_2D<Patch2> const& spline_patch_2)
{
    double max_error = 0;
    Kokkos::parallel_reduce(
            eval_coords.extent(0),
            KOKKOS_LAMBDA(int i, double& err) {
                double const eval_function = evaluators.deriv_1_and_2(eval_coords(i), splines);
                double expected_function;
                if (Coord<R>(eval_coords(i)) < ddc::discrete_space<BSplinesR<1>>().rmax()) {
                    expected_function
                            = single_evaluator_1.deriv_1_and_2(eval_coords(i), spline_patch_1);
                } else {
                    expected_function
                            = single_evaluator_2.deriv_1_and_2(eval_coords(i), spline_patch_2);
                }
                err = Kokkos::max(abs(eval_function - expected_function), err);
            },
            Kokkos::Max<double>(max_error));
    EXPECT_LE(max_error, 1e-15);
};

// Test deriv() ..................................................................................
template <class InterestDim>
void test_deriv(
        DeviceMultipatchSplineRThetaEvaluator const& evaluators,
        SplineRThetaEvaluator<1, DeviceExecSpace> const& single_evaluator_1,
        SplineRThetaEvaluator<2, DeviceExecSpace> const& single_evaluator_2,
        Kokkos::View<Coord<R, Theta>*> const& eval_coords,
        MultipatchField<ConstSplineCoeffOnPatch_2D, Patch1, Patch2> const& splines,
        ConstSplineCoeffOnPatch_2D<Patch1> const& spline_patch_1,
        ConstSplineCoeffOnPatch_2D<Patch2> const& spline_patch_2)
{
    double max_error = 0;
    Kokkos::parallel_reduce(
            eval_coords.extent(0),
            KOKKOS_LAMBDA(int i, double& err) {
                double const eval_function = evaluators.deriv<InterestDim>(eval_coords(i), splines);
                double expected_function;
                if (Coord<R>(eval_coords(i)) < ddc::discrete_space<BSplinesR<1>>().rmax()) {
                    expected_function
                            = single_evaluator_1.deriv<InterestDim>(eval_coords(i), spline_patch_1);
                } else {
                    expected_function
                            = single_evaluator_2.deriv<InterestDim>(eval_coords(i), spline_patch_2);
                }
                err = Kokkos::max(abs(eval_function - expected_function), err);
            },
            Kokkos::Max<double>(max_error));
    EXPECT_LE(max_error, 1e-15);
};

} // namespace


/* -----------------------------------------------------------------------------------------------
    Test operator() for a single coordinate on host.
    --------------------------------------------------------------------------------------------*/
TEST_F(MultipatchSplineEvaluatorTest, HostEvaluateOnSingleCoord)
{
    // Definition of MultipatchSplineEvaluator2D
    PatchLocator<HostExecSpace> const
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);
    ConstantExtrapolationRuleOnion<PatchLocator<HostExecSpace>> extrapolation_rule(r1_min, r2_max);
    HostMultipatchSplineRThetaEvaluator const evaluators(patch_locator, extrapolation_rule);

    // Compare the evaluated functions with the expected functions.
    // Evaluation points
    // --- patch 1
    typename Patch1::Coord12 const eval_coord_1_on_patch_1(0.5, 0.0);
    typename Patch1::Coord12 const eval_coord_1_on_patch_2(1.1, M_PI);

    // --- patch 2
    typename Patch2::Coord12 const eval_coord_2_on_patch_1(0.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 const eval_coord_2_on_patch_2(1.9, 4 * M_PI);
    typename Patch2::Coord12 const eval_coord_2_outside(2.3, -M_PI);

    auto const function_1_coef_host = ddc::create_mirror_and_copy(function_1_coef);
    auto const function_2_coef_host = ddc::create_mirror_and_copy(function_2_coef);
    MultipatchField<ConstSplineCoeffOnPatch_2D_host, Patch1, Patch2> const
            splines_host(get_field(function_1_coef_host), get_field(function_2_coef_host));

    host_t<DConstField<IdxRange<BSplinesR<1>, BSplinesTheta<1>>>> const const_function_1_coef_host(
            get_const_field(function_1_coef_host));
    host_t<DConstField<IdxRange<BSplinesR<2>, BSplinesTheta<2>>>> const const_function_2_coef_host(
            get_const_field(function_2_coef_host));

    // Check assignement operator on host.
    double eval_function = evaluators(eval_coord_1_on_patch_1, splines_host);
    double expected_function
            = evaluator_1_host(eval_coord_1_on_patch_1, const_function_1_coef_host);
    EXPECT_NEAR(eval_function, expected_function, 1e-15);

    eval_function = evaluators(eval_coord_1_on_patch_2, splines_host);
    expected_function = evaluator_2_host(eval_coord_1_on_patch_2, const_function_2_coef_host);
    EXPECT_NEAR(eval_function, expected_function, 1e-15);

    eval_function = evaluators(eval_coord_2_on_patch_1, splines_host);
    expected_function = evaluator_1_host(eval_coord_2_on_patch_1, const_function_1_coef_host);
    EXPECT_NEAR(eval_function, expected_function, 1e-15);

    eval_function = evaluators(eval_coord_2_on_patch_2, splines_host);
    expected_function = evaluator_2_host(eval_coord_2_on_patch_2, const_function_2_coef_host);
    EXPECT_NEAR(eval_function, expected_function, 1e-15);

    eval_function = evaluators(eval_coord_2_outside, splines_host);
    expected_function = evaluator_2_host(eval_coord_2_outside, const_function_2_coef_host);
    EXPECT_NEAR(eval_function, expected_function, 1e-15);

    /*
        REMARK: here Patch1::Coord12 == Patch2::Coord12 == Coord<R, Theta>.
            So we do not need conversion to apply the spline evaluators of the
            patch 1 a coordinate from the patch 2, and vice versa.
    */
}


/* -----------------------------------------------------------------------------------------------
    Test operator() for a single coordinate on device.
    --------------------------------------------------------------------------------------------*/
TEST_F(MultipatchSplineEvaluatorTest, DeviceEvaluateOnSingleCoord)
{
    // Definition of MultipatchSplineEvaluator2D
    PatchLocator<DeviceExecSpace> const
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);
    ConstantExtrapolationRuleOnion<PatchLocator<DeviceExecSpace>>
            extrapolation_rule(r1_min, r2_max);
    DeviceMultipatchSplineRThetaEvaluator const evaluators(patch_locator, extrapolation_rule);

    // Compare the evaluated functions with the expected functions.
    // Evaluation points
    // --- patch 1
    typename Patch1::Coord12 const eval_coord_1_on_patch_1(0.5, 0.0);
    typename Patch1::Coord12 const eval_coord_1_on_patch_2(1.1, M_PI);

    // --- patch 2
    typename Patch2::Coord12 const eval_coord_2_on_patch_1(0.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 const eval_coord_2_on_patch_2(1.9, 4 * M_PI);
    typename Patch2::Coord12 const eval_coord_2_outside(2.3, -M_PI);

    // Coordinates and patch indices on device.
    Kokkos::View<Coord<R, Theta>*, Kokkos::DefaultExecutionSpace::memory_space> coords("coords", 5);
    Kokkos::View<Coord<R, Theta>*, Kokkos::HostSpace> coords_host("coords_host", 5);

    coords_host(0) = eval_coord_1_on_patch_1;
    coords_host(1) = eval_coord_1_on_patch_2;
    coords_host(2) = eval_coord_2_on_patch_1;
    coords_host(3) = eval_coord_2_on_patch_2;
    coords_host(4) = eval_coord_2_outside;

    Kokkos::deep_copy(coords, coords_host);

    // Check assignement operator on device.
    test_operator_assignement(
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
TEST_F(MultipatchSplineEvaluatorTest, EvaluateOnCoordField)
{
    // Evaluation points
    // --- patch 1
    Patch1::IdxRange1 const reduced_idx_range_r1(
            Patch1::IdxRange1(Patch1::Idx1(0), Patch1::IdxStep1(idx_range_r1.size() - 1)));
    Patch1::IdxRange12 const reduced_idx_range_rtheta1(reduced_idx_range_r1, idx_range_theta1);
    FieldMem<Patch1::Coord12, Patch1::IdxRange12> eval_points_1_alloc(reduced_idx_range_rtheta1);
    Field<Patch1::Coord12, Patch1::IdxRange12> eval_points_1 = get_field(eval_points_1_alloc);

    // --- patch 1
    Patch2::IdxRange1 const reduced_idx_range_r2(
            Patch2::IdxRange1(Patch2::Idx1(0), Patch2::IdxStep1(idx_range_r2.size() - 1)));
    Patch2::IdxRange12 const reduced_idx_range_rtheta2(reduced_idx_range_r2, idx_range_theta2);
    FieldMem<Patch2::Coord12, Patch2::IdxRange12> eval_points_2_alloc(reduced_idx_range_rtheta2);
    Field<Patch2::Coord12, Patch2::IdxRange12> eval_points_2 = get_field(eval_points_2_alloc);

    set_eval_points_2D<Patch1::Grid1, Patch1::Grid2>(eval_points_1);
    set_eval_points_2D<Patch2::Grid1, Patch2::Grid2>(eval_points_2);

    // --- collection
    MultipatchField<CoordConstFieldOnPatch, Patch1, Patch2> const
            eval_points(get_const_field(eval_points_1), get_const_field(eval_points_2));


    // Evaluated functions
    // --- patch 1
    DFieldMem<Patch1::IdxRange12> eval_function_1_alloc(reduced_idx_range_rtheta1);
    DField<Patch1::IdxRange12> eval_function_1 = get_field(eval_function_1_alloc);

    // --- patch 2
    DFieldMem<Patch2::IdxRange12> eval_function_2_alloc(reduced_idx_range_rtheta2);
    DField<Patch2::IdxRange12> eval_function_2 = get_field(eval_function_2_alloc);

    // --- collection
    MultipatchField<DFieldOnPatch, Patch1, Patch2> const
            eval_functions(eval_function_1, eval_function_2);


    // Definition of MultipatchSplineEvaluator2D
    PatchLocator<DeviceExecSpace> const
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);
    ConstantExtrapolationRuleOnion<PatchLocator<DeviceExecSpace>>
            extrapolation_rule(r1_min, r2_max);
    DeviceMultipatchSplineRThetaEvaluator const evaluators(patch_locator, extrapolation_rule);

    // Evaluate the functions at the evaluation points.
    evaluators(eval_functions, eval_points, splines);

    // Expected functions
    // --- patch 1
    DFieldMem<Patch1::IdxRange12> expected_function_1_alloc(reduced_idx_range_rtheta1);
    DField<Patch1::IdxRange12> expected_function_1 = get_field(expected_function_1_alloc);

    // --- patch 2
    DFieldMem<Patch2::IdxRange12> expected_function_2_alloc(reduced_idx_range_rtheta2);
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
    ddc::for_each(reduced_idx_range_rtheta1, [&](typename Patch1::Idx12 const idx) {
        EXPECT_NEAR(eval_function_1_host(idx), expected_function_1_host(idx), 1e-15);
    });
    ddc::for_each(reduced_idx_range_rtheta2, [&](typename Patch2::Idx12 const idx) {
        EXPECT_NEAR(eval_function_2_host(idx), expected_function_2_host(idx), 1e-15);
    });
}


/* -----------------------------------------------------------------------------------------------
    Test deriv_dim_1(), deriv_dim_2(), deriv_1_and_2() and deriv<InterestDim>()
    for a single coordinate on host.
    --------------------------------------------------------------------------------------------*/
TEST_F(MultipatchSplineEvaluatorTest, HostDerivativesOnSingleCoord)
{
    // Definition of MultipatchSplineEvaluator2D
    PatchLocator<HostExecSpace> const
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);
    ConstantExtrapolationRuleOnion<PatchLocator<HostExecSpace>> extrapolation_rule(r1_min, r2_max);
    HostMultipatchSplineRThetaEvaluator const evaluators(patch_locator, extrapolation_rule);

    // Compare the evaluated functions with the expected functions.
    // Evaluation points
    // --- patch 1
    typename Patch1::Coord12 const eval_coord_1_on_patch_1(0.5, 0.0);
    typename Patch1::Coord12 const eval_coord_1_on_patch_2(1.1, M_PI);

    // --- patch 2
    typename Patch2::Coord12 const eval_coord_2_on_patch_1(0.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 const eval_coord_2_on_patch_2(1.9, 0.0);


    auto const function_1_coef_host = ddc::create_mirror_and_copy(function_1_coef);
    auto const function_2_coef_host = ddc::create_mirror_and_copy(function_2_coef);
    MultipatchField<ConstSplineCoeffOnPatch_2D_host, Patch1, Patch2> const
            splines_host(get_field(function_1_coef_host), get_field(function_2_coef_host));

    host_t<DConstField<IdxRange<BSplinesR<1>, BSplinesTheta<1>>>> const cst_function_1_coef_h(
            get_const_field(function_1_coef_host));
    host_t<DConstField<IdxRange<BSplinesR<2>, BSplinesTheta<2>>>> const cst_function_2_coef_h(
            get_const_field(function_2_coef_host));


    // Check deriv_dim_1(), deriv_dim_2(), deriv_1_and_2(), deriv() on host.
    // --- derivative 1
    double eval_deriv = evaluators.deriv_dim_1(eval_coord_1_on_patch_1, splines_host);
    double expected_deriv
            = evaluator_1_host.deriv_dim_1(eval_coord_1_on_patch_1, cst_function_1_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.deriv_dim_1(eval_coord_1_on_patch_2, splines_host);
    expected_deriv = evaluator_2_host.deriv_dim_1(eval_coord_1_on_patch_2, cst_function_2_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.deriv_dim_1(eval_coord_2_on_patch_1, splines_host);
    expected_deriv = evaluator_1_host.deriv_dim_1(eval_coord_2_on_patch_1, cst_function_1_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.deriv_dim_1(eval_coord_2_on_patch_2, splines_host);
    expected_deriv = evaluator_2_host.deriv_dim_1(eval_coord_2_on_patch_2, cst_function_2_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);


    // --- derivative 2
    eval_deriv = evaluators.deriv_dim_2(eval_coord_1_on_patch_1, splines_host);
    expected_deriv = evaluator_1_host.deriv_dim_2(eval_coord_1_on_patch_1, cst_function_1_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.deriv_dim_2(eval_coord_1_on_patch_2, splines_host);
    expected_deriv = evaluator_2_host.deriv_dim_2(eval_coord_1_on_patch_2, cst_function_2_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.deriv_dim_2(eval_coord_2_on_patch_1, splines_host);
    expected_deriv = evaluator_1_host.deriv_dim_2(eval_coord_2_on_patch_1, cst_function_1_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.deriv_dim_2(eval_coord_2_on_patch_2, splines_host);
    expected_deriv = evaluator_2_host.deriv_dim_2(eval_coord_2_on_patch_2, cst_function_2_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);


    // --- cross-derivative
    eval_deriv = evaluators.deriv_1_and_2(eval_coord_1_on_patch_1, splines_host);
    expected_deriv = evaluator_1_host.deriv_1_and_2(eval_coord_1_on_patch_1, cst_function_1_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.deriv_1_and_2(eval_coord_1_on_patch_2, splines_host);
    expected_deriv = evaluator_2_host.deriv_1_and_2(eval_coord_1_on_patch_2, cst_function_2_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.deriv_1_and_2(eval_coord_2_on_patch_1, splines_host);
    expected_deriv = evaluator_1_host.deriv_1_and_2(eval_coord_2_on_patch_1, cst_function_1_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.deriv_1_and_2(eval_coord_2_on_patch_2, splines_host);
    expected_deriv = evaluator_2_host.deriv_1_and_2(eval_coord_2_on_patch_2, cst_function_2_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);


    // --- template derivative
    eval_deriv = evaluators.deriv<R>(eval_coord_1_on_patch_1, splines_host);
    expected_deriv = evaluator_1_host.deriv<R>(eval_coord_1_on_patch_1, cst_function_1_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.deriv<R>(eval_coord_1_on_patch_2, splines_host);
    expected_deriv = evaluator_2_host.deriv<R>(eval_coord_1_on_patch_2, cst_function_2_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.template deriv<Theta>(eval_coord_2_on_patch_1, splines_host);
    expected_deriv = evaluator_1_host.deriv<Theta>(eval_coord_2_on_patch_1, cst_function_1_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);

    eval_deriv = evaluators.template deriv<Theta>(eval_coord_2_on_patch_2, splines_host);
    expected_deriv = evaluator_2_host.deriv<Theta>(eval_coord_2_on_patch_2, cst_function_2_coef_h);
    EXPECT_NEAR(eval_deriv, expected_deriv, 1e-15);
}


/* -----------------------------------------------------------------------------------------------
    Test deriv_dim_1(), deriv_dim_2(), deriv_1_and_2() and deriv<InterestDim>()
    for a single coordinate on device.
    --------------------------------------------------------------------------------------------*/
TEST_F(MultipatchSplineEvaluatorTest, DeviceDerivativesOnSingleCoord)
{
    // Definition of MultipatchSplineEvaluator2D
    PatchLocator<DeviceExecSpace> const
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);
    ConstantExtrapolationRuleOnion<PatchLocator<DeviceExecSpace>>
            extrapolation_rule(r1_min, r2_max);
    DeviceMultipatchSplineRThetaEvaluator const evaluators(patch_locator, extrapolation_rule);

    // Compare the evaluated functions with the expected functions.
    // Evaluation points
    // --- patch 1
    typename Patch1::Coord12 const eval_coord_1_on_patch_1(0.5, 0.0);
    typename Patch1::Coord12 const eval_coord_1_on_patch_2(1.1, M_PI);

    // --- patch 2
    typename Patch2::Coord12 const eval_coord_2_on_patch_1(0.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 const eval_coord_2_on_patch_2(1.9, 0.0);

    Kokkos::View<Coord<R, Theta>*> coords("coords", 4);
    Kokkos::View<Coord<R, Theta>*, Kokkos::HostSpace> coords_host("coords_host", 4);

    coords_host(0) = eval_coord_1_on_patch_1;
    coords_host(1) = eval_coord_1_on_patch_2;
    coords_host(2) = eval_coord_2_on_patch_1;
    coords_host(3) = eval_coord_2_on_patch_2;

    Kokkos::deep_copy(coords, coords_host);

    // Check deriv_dim_1(), deriv_dim_2(), deriv_1_and_2(), deriv() on device.
    // --- derivativative 1 ----------------------
    test_deriv_dim_1(
            evaluators,
            evaluator_1,
            evaluator_2,
            coords,
            splines,
            get_const_field(function_1_coef),
            get_const_field(function_2_coef));

    // --- derivativative 2 ----------------------
    test_deriv_dim_2(
            evaluators,
            evaluator_1,
            evaluator_2,
            coords,
            splines,
            get_const_field(function_1_coef),
            get_const_field(function_2_coef));

    // --- cross-derivativative ------------------
    test_deriv_1_and_2(
            evaluators,
            evaluator_1,
            evaluator_2,
            coords,
            splines,
            get_const_field(function_1_coef),
            get_const_field(function_2_coef));

    // --- derivatives on template dimension -----
    test_deriv<R>(
            evaluators,
            evaluator_1,
            evaluator_2,
            coords,
            splines,
            get_const_field(function_1_coef),
            get_const_field(function_2_coef));

    test_deriv<Theta>(
            evaluators,
            evaluator_1,
            evaluator_2,
            coords,
            splines,
            get_const_field(function_1_coef),
            get_const_field(function_2_coef));
}


/* -----------------------------------------------------------------------------------------------
    Test deriv_dim_1(), deriv_dim_2(), deriv_1_and_2() and deriv<InterestDim>()
    for single coordinates outside of the domain. The operators are supposed to
    return an error.
    --------------------------------------------------------------------------------------------*/
TEST_F(MultipatchSplineEvaluatorTest, DerivativesOnSingleCoordDeathTest)
{
    // Definition of MultipatchSplineEvaluator2D
    PatchLocator<DeviceExecSpace> const
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);
    ConstantExtrapolationRuleOnion<PatchLocator<DeviceExecSpace>>
            extrapolation_rule(r1_min, r2_max);
    DeviceMultipatchSplineRThetaEvaluator const evaluators(patch_locator, extrapolation_rule);

    // Evaluation points
    typename Patch1::Coord12 const eval_coord_1_outside(2.5, 0.0);
    typename Patch2::Coord12 const eval_coord_2_outside(2.5, 0.0);

    // Test the Kokkos::abort only on host.
    if constexpr (std::is_same_v<
                          Kokkos::DefaultExecutionSpace,
                          Kokkos::DefaultHostExecutionSpace>) {
        // --- derivativative 1
        EXPECT_DEATH(
                evaluators.deriv_dim_1(eval_coord_1_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");

        EXPECT_DEATH(
                evaluators.deriv_dim_1(eval_coord_2_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");

        // --- derivativative 2
        EXPECT_DEATH(
                evaluators.deriv_dim_2(eval_coord_1_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");

        EXPECT_DEATH(
                evaluators.deriv_dim_2(eval_coord_2_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");

        // --- cross-derivativative
        EXPECT_DEATH(
                evaluators.deriv_1_and_2(eval_coord_1_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");

        EXPECT_DEATH(
                evaluators.deriv_1_and_2(eval_coord_2_outside, splines),
                "The evaluation coordinate has to be on a patch."
                "No extrapolation rule for derivatives. \n");
    }
}


/* -----------------------------------------------------------------------------------------------
    Test deriv_dim_1(), deriv_dim_2(), deriv_1_and_2() and deriv<InterestDim>()
    for fields of coordinates called from host and stored on device.
    --------------------------------------------------------------------------------------------*/
TEST_F(MultipatchSplineEvaluatorTest, DerivativativesOnCoordField)
{
    // Evaluation points
    // --- patch 1
    Patch1::IdxRange1 const reduced_idx_range_r1(
            Patch1::IdxRange1(Patch1::Idx1(0), Patch1::IdxStep1(idx_range_r1.size() - 1)));
    Patch1::IdxRange2 const reduced_idx_range_theta1(
            Patch1::IdxRange2(Patch1::Idx2(0), Patch1::IdxStep2(idx_range_theta1.size() - 1)));
    Patch1::IdxRange12 const
            reduced_idx_range_rtheta1(reduced_idx_range_r1, reduced_idx_range_theta1);
    FieldMem<Patch1::Coord12, Patch1::IdxRange12> eval_points_1_alloc(reduced_idx_range_rtheta1);
    Field<Patch1::Coord12, Patch1::IdxRange12> eval_points_1 = get_field(eval_points_1_alloc);

    // --- patch 1
    Patch2::IdxRange1 const reduced_idx_range_r2(
            Patch2::IdxRange1(Patch2::Idx1(0), Patch2::IdxStep1(idx_range_r2.size() - 1)));
    Patch2::IdxRange2 const reduced_idx_range_theta2(
            Patch2::IdxRange2(Patch2::Idx2(0), Patch2::IdxStep2(idx_range_theta2.size() - 1)));
    Patch2::IdxRange12 const
            reduced_idx_range_rtheta2(reduced_idx_range_r2, reduced_idx_range_theta2);
    FieldMem<Patch2::Coord12, Patch2::IdxRange12> eval_points_2_alloc(reduced_idx_range_rtheta2);
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
    DFieldMem<Patch1::IdxRange12> eval_derivs_1_patch_1_alloc(reduced_idx_range_rtheta1);
    DField<Patch1::IdxRange12> eval_derivs_1_patch_1 = get_field(eval_derivs_1_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> eval_derivs_1_patch_2_alloc(reduced_idx_range_rtheta2);
    DField<Patch2::IdxRange12> eval_derivs_1_patch_2 = get_field(eval_derivs_1_patch_2_alloc);

    // --- derivativative 2
    DFieldMem<Patch1::IdxRange12> eval_derivs_2_patch_1_alloc(reduced_idx_range_rtheta1);
    DField<Patch1::IdxRange12> eval_derivs_2_patch_1 = get_field(eval_derivs_2_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> eval_derivs_2_patch_2_alloc(reduced_idx_range_rtheta2);
    DField<Patch2::IdxRange12> eval_derivs_2_patch_2 = get_field(eval_derivs_2_patch_2_alloc);

    // --- cross-derivatives
    DFieldMem<Patch1::IdxRange12> eval_derivs_12_patch_1_alloc(reduced_idx_range_rtheta1);
    DField<Patch1::IdxRange12> eval_derivs_12_patch_1 = get_field(eval_derivs_12_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> eval_derivs_12_patch_2_alloc(reduced_idx_range_rtheta2);
    DField<Patch2::IdxRange12> eval_derivs_12_patch_2 = get_field(eval_derivs_12_patch_2_alloc);

    // ----- collections
    MultipatchField<DFieldOnPatch, Patch1, Patch2> const
            eval_derivs_1(eval_derivs_1_patch_1, eval_derivs_1_patch_2);
    MultipatchField<DFieldOnPatch, Patch1, Patch2> const
            eval_derivs_2(eval_derivs_2_patch_1, eval_derivs_2_patch_2);
    MultipatchField<DFieldOnPatch, Patch1, Patch2> const
            eval_derivs_12(eval_derivs_12_patch_1, eval_derivs_12_patch_2);


    // Definition of MultipatchSplineEvaluator2D
    PatchLocator<DeviceExecSpace> const
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);
    ConstantExtrapolationRuleOnion<PatchLocator<DeviceExecSpace>>
            extrapolation_rule(r1_min, r2_max);
    DeviceMultipatchSplineRThetaEvaluator const evaluators(patch_locator, extrapolation_rule);

    // Evaluate the functions at the evaluation points.
    evaluators.deriv_dim_1(eval_derivs_1, eval_points, splines);
    evaluators.deriv_dim_2(eval_derivs_2, eval_points, splines);
    evaluators.deriv_1_and_2(eval_derivs_12, eval_points, splines);


    // Expected functions
    // --- derivative 1
    DFieldMem<Patch1::IdxRange12> expected_derivs_1_patch_1_alloc(reduced_idx_range_rtheta1);
    DField<Patch1::IdxRange12> expected_derivs_1_patch_1
            = get_field(expected_derivs_1_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> expected_derivs_1_patch_2_alloc(reduced_idx_range_rtheta2);
    DField<Patch2::IdxRange12> expected_derivs_1_patch_2
            = get_field(expected_derivs_1_patch_2_alloc);

    // --- derivative 2
    DFieldMem<Patch1::IdxRange12> expected_derivs_2_patch_1_alloc(reduced_idx_range_rtheta1);
    DField<Patch1::IdxRange12> expected_derivs_2_patch_1
            = get_field(expected_derivs_2_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> expected_derivs_2_patch_2_alloc(reduced_idx_range_rtheta2);
    DField<Patch2::IdxRange12> expected_derivs_2_patch_2
            = get_field(expected_derivs_2_patch_2_alloc);

    // --- cross-derivative
    DFieldMem<Patch1::IdxRange12> expected_derivs_12_patch_1_alloc(reduced_idx_range_rtheta1);
    DField<Patch1::IdxRange12> expected_derivs_12_patch_1
            = get_field(expected_derivs_12_patch_1_alloc);

    DFieldMem<Patch2::IdxRange12> expected_derivs_12_patch_2_alloc(reduced_idx_range_rtheta2);
    DField<Patch2::IdxRange12> expected_derivs_12_patch_2
            = get_field(expected_derivs_12_patch_2_alloc);


    DConstField<typename Patch1::IdxRangeBS12> const const_function_1_coef
            = get_const_field(function_1_coef);
    DConstField<typename Patch2::IdxRangeBS12> const const_function_2_coef
            = get_const_field(function_2_coef);


    evaluator_1.deriv_dim_1(expected_derivs_1_patch_1, const_eval_points_1, const_function_1_coef);
    evaluator_2.deriv_dim_1(expected_derivs_1_patch_2, const_eval_points_2, const_function_2_coef);

    evaluator_1.deriv_dim_2(expected_derivs_2_patch_1, const_eval_points_1, const_function_1_coef);
    evaluator_2.deriv_dim_2(expected_derivs_2_patch_2, const_eval_points_2, const_function_2_coef);

    evaluator_1
            .deriv_1_and_2(expected_derivs_12_patch_1, const_eval_points_1, const_function_1_coef);
    evaluator_2
            .deriv_1_and_2(expected_derivs_12_patch_2, const_eval_points_2, const_function_2_coef);


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
    ddc::for_each(reduced_idx_range_rtheta1, [&](typename Patch1::Idx12 const idx) {
        EXPECT_NEAR(expected_derivs_1_patch_1_host(idx), eval_derivs_1_patch_1_host(idx), 1e-15);
        EXPECT_NEAR(expected_derivs_2_patch_1_host(idx), eval_derivs_2_patch_1_host(idx), 1e-15);
        EXPECT_NEAR(expected_derivs_12_patch_1_host(idx), eval_derivs_12_patch_1_host(idx), 1e-15);
    });
    ddc::for_each(reduced_idx_range_rtheta2, [&](typename Patch2::Idx12 const idx) {
        EXPECT_NEAR(expected_derivs_1_patch_2_host(idx), eval_derivs_1_patch_2_host(idx), 1e-15);
        EXPECT_NEAR(expected_derivs_2_patch_2_host(idx), eval_derivs_2_patch_2_host(idx), 1e-15);
        EXPECT_NEAR(expected_derivs_12_patch_2_host(idx), eval_derivs_12_patch_2_host(idx), 1e-15);
    });
}


/* -----------------------------------------------------------------------------------------------
    Test integrate() on host.
    --------------------------------------------------------------------------------------------*/
TEST_F(MultipatchSplineEvaluatorTest, HostIntegrateOnCoordField)
{
    auto function_1_coef_host = ddc::create_mirror_and_copy(function_1_coef);
    auto function_2_coef_host = ddc::create_mirror_and_copy(function_2_coef);
    MultipatchField<ConstSplineCoeffOnPatch_2D_host, Patch1, Patch2> const splines_host(
            get_const_field(function_1_coef_host),
            get_const_field(function_2_coef_host));

    // Definition of MultipatchSplineEvaluator2D
    PatchLocator<HostExecSpace> const
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);
    ConstantExtrapolationRuleOnion<PatchLocator<HostExecSpace>> extrapolation_rule(r1_min, r2_max);
    HostMultipatchSplineRThetaEvaluator const evaluators(patch_locator, extrapolation_rule);

    // Evaluate the functions at the evaluation points.
    Kokkos::View<double[2], HostExecSpace> eval_integrals("eval_integrals", 2);
    evaluators.integrate(eval_integrals, splines_host);

    // Expected functions
    double const expected_integral_1 = 0.0;
    double const expected_integral_2 = 0.0;

    // Compare the evaluated functions with the expected functions.
    EXPECT_NEAR(eval_integrals(0), expected_integral_1, 1e-15);
    EXPECT_NEAR(eval_integrals(1), expected_integral_2, 1e-15);

    /*
        REMARK: Cannot compare to SplineEvaluator2D as long as the .integrate()
            function is only defined on batch domain.
    */
}
