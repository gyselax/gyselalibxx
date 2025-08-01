
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "multipatch_field.hpp"
#include "types.hpp"


/// @brief Convert a local coordinate into a global coordinate.
template <typename Xg, typename Yg, typename X_loc, typename Y_loc>
Coord<Xg, Yg> get_global_coord(Coord<X_loc, Y_loc> const& local_coord)
{
    double const x = ddc::select<X_loc>(local_coord);
    double const y = ddc::select<Y_loc>(local_coord);
    return Coord<Xg, Yg>(x, y);
}


/// @brief Initialise the function with f(x,y) = cos(2/3*pi*x)sin(y).
template <class Grid1, class Grid2, class Layout>
void initialise_2D_function(DField<IdxRange<Grid1, Grid2>, Kokkos::HostSpace, Layout> function)
{
    ddc::for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        double const xg = ddc::coordinate(Idx<Grid1>(idx));
        double const yg = ddc::coordinate(Idx<Grid2>(idx));
        function(idx) = Kokkos::cos(xg * 2. / 3. * M_PI) * Kokkos::sin(yg);
    });
}

/// @brief Initialise all the functions defined on the patches.
template <class... Patches>
void initialise_all_functions(
        MultipatchField<DerivFieldOnPatch_host, Patches...> const& functions_and_derivs)
{
    (initialise_2D_function<typename Patches::Grid1, typename Patches::Grid2>(
             (functions_and_derivs.template get<Patches>()).get_values_field()),
     ...);
}



// template <
//         class Patch,
//         ddc::BoundCond BoundCondXmin,
//         ddc::BoundCond BoundCondXmax,
//         ddc::BoundCond BoundCondYmin,
//         ddc::BoundCond BoundCondYmax,
//         class BSplinesXg,
//         class BSplinesYg,
//         class SplineRThetagEvaluator>
// void check_local_and_global_splines_agreement(
//         typename Patch::IdxRange12 const& idx_range_xy,
//         IdxRange1SliceOnPatch<Patch> const& idx_range_slice_x,
//         IdxRange2SliceOnPatch<Patch> const& idx_range_slice_y,
//         DerivFieldOnPatch_host<Patch> function_and_derivs,
//         SplineRThetagEvaluator const& evaluator_g,
//         host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
// {
//     using HostExecSpace = Kokkos::DefaultHostExecutionSpace;

//     using DimX = typename Patch::Dim1;
//     using DimY = typename Patch::Dim2;
//     using DerivX = typename ddc::Deriv<DimX>;
//     using DerivY = typename ddc::Deriv<DimY>;
//     using GridX = typename Patch::Grid1;
//     using GridY = typename Patch::Grid2;

//     using Xg = typename BSplinesXg::continuous_dimension_type;
//     using Yg = typename BSplinesYg::continuous_dimension_type;

//     // Build the local spline representation -------------------------------------------------
//     ddc::SplineBuilder2D<
//             HostExecSpace,
//             typename HostExecSpace::memory_space,
//             typename Patch::BSplines1,
//             typename Patch::BSplines2,
//             typename Patch::Grid1,
//             typename Patch::Grid2,
//             BoundCondXmin,
//             BoundCondXmax,
//             BoundCondYmin,
//             BoundCondYmax,
//             ddc::SplineSolver::LAPACK>
//             builder(idx_range_xy);

//     SplineCoeffMemOnPatch_2D_host<Patch> function_coef_alloc(
//             builder.batched_spline_domain(idx_range_xy));
//     SplineCoeffOnPatch_2D_host<Patch> function_coef = get_field(function_coef_alloc);

//     // Get the fields on the right layout.
//     // --- extract each fields.
//     // ------ function.
//     DField<IdxRange<GridX, GridY>, Kokkos::HostSpace, Kokkos::layout_stride> function_extracted
//             = function_and_derivs.get_values_field();

//     // ------ first derivatives.
//     IdxRange<DerivX> idx_range_deriv_x(Idx<DerivX>(1), IdxStep<DerivX>(1));
//     IdxRange<DerivY> idx_range_deriv_y(Idx<DerivY>(1), IdxStep<DerivY>(1));
//     IdxRange<DerivX, DerivY> idx_range_deriv_x_deriv_y(idx_range_deriv_x, idx_range_deriv_y);

//     IdxRange<GridX> idx_range_x(idx_range_xy);
//     IdxRange<GridY> idx_range_y(idx_range_xy);

//     Idx<GridX> idx_slice_xmin(idx_range_slice_x.front());
//     Idx<GridX> idx_slice_xmax(idx_range_slice_x.back());
//     Idx<GridY> idx_slice_ymin(idx_range_slice_y.front());
//     Idx<GridY> idx_slice_ymax(idx_range_slice_y.back());

//     IdxRange<GridX> idx_range_slice_xmin(idx_slice_xmin, IdxStep<GridX>(1));
//     IdxRange<GridX> idx_range_slice_xmax(idx_slice_xmax, IdxStep<GridX>(1));
//     IdxRange<GridY> idx_range_slice_ymin(idx_slice_ymin, IdxStep<GridY>(1));
//     IdxRange<GridY> idx_range_slice_ymax(idx_slice_ymax, IdxStep<GridY>(1));

//     IdxRange<DerivX, GridX, GridY>
//             idx_range_deriv_xmin(idx_range_deriv_x, idx_range_slice_xmin, idx_range_y);
//     IdxRange<DerivX, GridX, GridY>
//             idx_range_deriv_xmax(idx_range_deriv_x, idx_range_slice_xmax, idx_range_y);
//     IdxRange<GridX, DerivY, GridY>
//             idx_range_deriv_ymin(idx_range_x, idx_range_deriv_y, idx_range_slice_ymin);
//     IdxRange<GridX, DerivY, GridY>
//             idx_range_deriv_ymax(idx_range_x, idx_range_deriv_y, idx_range_slice_ymax);

//     DField<IdxRange<DerivX, GridY>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_xmin_extracted
//             = function_and_derivs[idx_range_deriv_xmin][idx_slice_xmin];
//     DField<IdxRange<DerivX, GridY>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_xmax_extracted
//             = function_and_derivs[idx_range_deriv_xmax][idx_slice_xmax];

//     DField<IdxRange<GridX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_ymin_extracted
//             = function_and_derivs[idx_range_deriv_ymin][idx_slice_ymin];
//     DField<IdxRange<GridX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_ymax_extracted
//             = function_and_derivs[idx_range_deriv_ymax][idx_slice_ymax];

//     // ------ cross-derivatives.
//     IdxRange<DerivX, GridX, DerivY, GridY> idx_range_deriv_xy_min_min(
//             idx_range_deriv_x,
//             idx_range_slice_xmin,
//             idx_range_deriv_y,
//             idx_range_slice_ymin);
//     IdxRange<DerivX, GridX, DerivY, GridY> idx_range_deriv_xy_max_min(
//             idx_range_deriv_x,
//             idx_range_slice_xmax,
//             idx_range_deriv_y,
//             idx_range_slice_ymin);
//     IdxRange<DerivX, GridX, DerivY, GridY> idx_range_deriv_xy_min_max(
//             idx_range_deriv_x,
//             idx_range_slice_xmin,
//             idx_range_deriv_y,
//             idx_range_slice_ymax);
//     IdxRange<DerivX, GridX, DerivY, GridY> idx_range_deriv_xy_max_max(
//             idx_range_deriv_x,
//             idx_range_slice_xmax,
//             idx_range_deriv_y,
//             idx_range_slice_ymax);

//     DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
//             deriv_xy_min_min_extracted
//             = function_and_derivs[idx_range_deriv_xy_min_min][idx_slice_xmin][idx_slice_ymin];
//     DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
//             deriv_xy_max_min_extracted
//             = function_and_derivs[idx_range_deriv_xy_max_min][idx_slice_xmax][idx_slice_ymin];
//     DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
//             deriv_xy_min_max_extracted
//             = function_and_derivs[idx_range_deriv_xy_min_max][idx_slice_xmin][idx_slice_ymax];
//     DField<IdxRange<DerivX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride>
//             deriv_xy_max_max_extracted
//             = function_and_derivs[idx_range_deriv_xy_max_max][idx_slice_xmax][idx_slice_ymax];


//     // --- define fields on the correct layout.
//     // ------ allocate memory
//     host_t<DFieldMem<IdxRange<GridX, GridY>>> function_alloc(idx_range_xy);

//     IdxRange<DerivX, GridY> idx_range_deriv_x_y(idx_range_deriv_x, idx_range_y);
//     host_t<DFieldMem<IdxRange<DerivX, GridY>>> deriv_xmin_alloc(idx_range_deriv_x_y);
//     host_t<DFieldMem<IdxRange<DerivX, GridY>>> deriv_xmax_alloc(idx_range_deriv_x_y);

//     IdxRange<GridX, DerivY> idx_range_x_deriv_y(idx_range_x, idx_range_deriv_y);
//     host_t<DFieldMem<IdxRange<GridX, DerivY>>> deriv_ymin_alloc(idx_range_x_deriv_y);
//     host_t<DFieldMem<IdxRange<GridX, DerivY>>> deriv_ymax_alloc(idx_range_x_deriv_y);

//     host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_min_min_alloc(idx_range_deriv_x_deriv_y);
//     host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_max_min_alloc(idx_range_deriv_x_deriv_y);
//     host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_min_max_alloc(idx_range_deriv_x_deriv_y);
//     host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_max_max_alloc(idx_range_deriv_x_deriv_y);

//     // ------ define span
//     host_t<DField<IdxRange<GridX, GridY>>> function(function_alloc);

//     host_t<DField<IdxRange<DerivX, GridY>>> deriv_xmin(deriv_xmin_alloc);
//     host_t<DField<IdxRange<DerivX, GridY>>> deriv_xmax(deriv_xmax_alloc);

//     host_t<DField<IdxRange<GridX, DerivY>>> deriv_ymin(deriv_ymin_alloc);
//     host_t<DField<IdxRange<GridX, DerivY>>> deriv_ymax(deriv_ymax_alloc);

//     host_t<DField<IdxRange<DerivX, DerivY>>> deriv_xy_min_min(deriv_xy_min_min_alloc);
//     host_t<DField<IdxRange<DerivX, DerivY>>> deriv_xy_max_min(deriv_xy_max_min_alloc);
//     host_t<DField<IdxRange<DerivX, DerivY>>> deriv_xy_min_max(deriv_xy_min_max_alloc);
//     host_t<DField<IdxRange<DerivX, DerivY>>> deriv_xy_max_max(deriv_xy_max_max_alloc);

//     // ------ initialise data from the fields on layout stride.
//     ddc::for_each(idx_range_xy, [&](Idx<GridX, GridY> const idx) {
//         function(idx) = function_extracted(idx);
//     });

//     ddc::for_each(idx_range_deriv_x_y, [&](Idx<DerivX, GridY> const idx) {
//         deriv_xmin(idx) = deriv_xmin_extracted(idx);
//         deriv_xmax(idx) = deriv_xmax_extracted(idx);
//     });

//     ddc::for_each(idx_range_x_deriv_y, [&](Idx<GridX, DerivY> const idx) {
//         deriv_ymin(idx) = deriv_ymin_extracted(idx);
//         deriv_ymax(idx) = deriv_ymax_extracted(idx);
//     });

//     Idx<DerivX, DerivY> const idx_cross_deriv(idx_range_deriv_x_deriv_y.front());
//     deriv_xy_min_min(idx_cross_deriv) = deriv_xy_min_min_extracted(idx_cross_deriv);
//     deriv_xy_max_min(idx_cross_deriv) = deriv_xy_max_min_extracted(idx_cross_deriv);
//     deriv_xy_min_max(idx_cross_deriv) = deriv_xy_min_max_extracted(idx_cross_deriv);
//     deriv_xy_max_max(idx_cross_deriv) = deriv_xy_max_max_extracted(idx_cross_deriv);


//     // If the boundary is not a ddc::BoundCond::HERMITE, we don't use derivatives.
//     if constexpr (
//             (BoundCondYmin == ddc::BoundCond::HERMITE)
//             && (BoundCondYmax == ddc::BoundCond::HERMITE)) {
//         builder(function_coef,
//                 get_const_field(function),
//                 std::optional(get_const_field(deriv_xmin)),
//                 std::optional(get_const_field(deriv_xmax)),
//                 std::optional(get_const_field(deriv_ymin)),
//                 std::optional(get_const_field(deriv_ymax)),
//                 std::optional(get_const_field(deriv_xy_min_min)),
//                 std::optional(get_const_field(deriv_xy_max_min)),
//                 std::optional(get_const_field(deriv_xy_min_max)),
//                 std::optional(get_const_field(deriv_xy_max_max)));
//     } else if constexpr (
//             (BoundCondYmin == ddc::BoundCond::HERMITE)
//             && !(BoundCondYmax == ddc::BoundCond::HERMITE)) {
//         builder(function_coef,
//                 get_const_field(function),
//                 std::optional(get_const_field(deriv_xmin)),
//                 std::optional(get_const_field(deriv_xmax)),
//                 std::optional(get_const_field(deriv_ymin)),
//                 std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
//                 std::optional(get_const_field(deriv_xy_min_min)),
//                 std::optional(get_const_field(deriv_xy_max_min)),
//                 std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
//                 std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt});
//     } else if constexpr (
//             !(BoundCondYmin == ddc::BoundCond::HERMITE)
//             && (BoundCondYmax == ddc::BoundCond::HERMITE)) {
//         builder(function_coef,
//                 get_const_field(function),
//                 std::optional(get_const_field(deriv_xmin)),
//                 std::optional(get_const_field(deriv_xmax)),
//                 std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
//                 std::optional(get_const_field(deriv_ymax)),
//                 std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
//                 std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
//                 std::optional(get_const_field(deriv_xy_min_max)),
//                 std::optional(get_const_field(deriv_xy_max_max)));
//     } else {
//         builder(function_coef,
//                 get_const_field(function),
//                 std::optional(get_const_field(deriv_xmin)),
//                 std::optional(get_const_field(deriv_xmax)),
//                 std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
//                 std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
//                 std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
//                 std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
//                 std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
//                 std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt});
//     }


//     // Define local spline evaluator ---------------------------------------------------------
//     Coord<DimX> const x_min(ddc::discrete_space<typename Patch::BSplines1>().rmin());
//     Coord<DimX> const x_max(ddc::discrete_space<typename Patch::BSplines1>().rmax());
//     Coord<DimY> const y_min(ddc::discrete_space<typename Patch::BSplines2>().rmin());
//     Coord<DimY> const y_max(ddc::discrete_space<typename Patch::BSplines2>().rmax());
//     ddc::ConstantExtrapolationRule<DimX, DimY> bc_xmin(x_min, y_min, y_max);
//     ddc::ConstantExtrapolationRule<DimX, DimY> bc_xmax(x_max, y_min, y_max);
//     ddc::ConstantExtrapolationRule<DimY, DimX> bc_ymin(y_min, x_min, x_max);
//     ddc::ConstantExtrapolationRule<DimY, DimX> bc_ymax(y_max, x_min, x_max);
//     ddc::SplineEvaluator2D<
//             HostExecSpace,
//             typename HostExecSpace::memory_space,
//             typename Patch::BSplines1,
//             typename Patch::BSplines2,
//             typename Patch::Grid1,
//             typename Patch::Grid2,
//             ddc::ConstantExtrapolationRule<DimX, DimY>,
//             ddc::ConstantExtrapolationRule<DimX, DimY>,
//             ddc::ConstantExtrapolationRule<DimY, DimX>,
//             ddc::ConstantExtrapolationRule<DimY, DimX>>
//             evaluator(bc_xmin, bc_xmax, bc_ymin, bc_ymax);

//     // Define evaluation points at the centre of the cells -----------------------------------
//     host_t<CoordFieldMemOnPatch<Patch>> eval_points_alloc(idx_range_xy);
//     host_t<CoordFieldOnPatch<Patch>> eval_points(eval_points_alloc);

//     // The evaluation points are placed in the middle of the cells, except for the last one.
//     ddc::for_each(idx_range_xy, [&](typename Patch::Idx12 const idx) {
//         Coord<DimX, DimY> const mesh_point(ddc::coordinate(idx));
//         typename Patch::Idx1 idx_1(idx);
//         typename Patch::Idx2 idx_2(idx);
//         Coord<DimX> dx = (idx_1 != typename Patch::IdxRange1(idx_range_xy).back())
//                          * distance_at_right(idx_1);
//         Coord<DimY> dy = (idx_2 != typename Patch::IdxRange2(idx_range_xy).back())
//                          * distance_at_right(idx_2);
//         eval_points(idx) = mesh_point + Coord<DimX, DimY>(dx, dy);
//     });

//     // Evaluate and compare the local and global spline representations ----------------------
//     ddc::for_each(idx_range_xy, [&](typename Patch::Idx12 const idx) {
//         Coord<Xg, Yg> const eval_point_g(get_global_coord<Xg,Yg>(eval_points(idx)));

//         double local_spline = evaluator(eval_points(idx), get_const_field(function_coef));
//         double global_spline = evaluator_g(eval_point_g, get_const_field(function_g_coef));

//         EXPECT_NEAR(local_spline, global_spline, 1e-14);
//     });
// }