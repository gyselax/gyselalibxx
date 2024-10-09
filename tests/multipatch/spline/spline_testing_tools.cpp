// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "spline_testing_tools.hpp"

DField<Patch1::IdxRangeBS12> MultipatchSplineOnionShapeTest::set_spline_1(
        DField<Patch1::IdxRangeBS12> const& function_1_coef)
{
    SplineRThetaBuilder<1, DeviceExecSpace> const builder_1(idx_range_rtheta1);

    // Exact function
    DFieldMem<Patch1::IdxRange12> function_1_alloc(idx_range_rtheta1);
    DField<Patch1::IdxRange12> function_1 = get_field(function_1_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_rtheta1,
            KOKKOS_LAMBDA(Patch1::Idx12 const idx) {
                double const r = ddc::coordinate(Patch1::Idx1(idx));
                double const theta = ddc::coordinate(Patch1::Idx2(idx));
                function_1(idx) = r * Kokkos::sin(theta);
            });

    // Build the spline representations on each patch
    builder_1(function_1_coef, get_const_field(function_1));
    return function_1_coef;
};


DField<Patch2::IdxRangeBS12> MultipatchSplineOnionShapeTest::set_spline_2(
        DField<Patch2::IdxRangeBS12> const& function_2_coef)
{
    SplineRThetaBuilder<2, DeviceExecSpace> const builder_2(idx_range_rtheta2);

    // Exact function
    DFieldMem<Patch2::IdxRange12> function_2_alloc(idx_range_rtheta2);
    DField<Patch2::IdxRange12> function_2 = get_field(function_2_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_rtheta2,
            KOKKOS_LAMBDA(Patch2::Idx12 const idx) {
                double const r = ddc::coordinate(Patch2::Idx1(idx));
                double const theta = ddc::coordinate(Patch2::Idx2(idx));
                function_2(idx) = r * r * Kokkos::sin(theta);
            });

    // Build the spline representations on each patch
    builder_2(function_2_coef, get_const_field(function_2));
    return function_2_coef;
};