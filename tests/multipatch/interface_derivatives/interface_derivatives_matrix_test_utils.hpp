
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "multipatch_field.hpp"
#include "types.hpp"

/*
    This file defined useful operators for the test on the InterfaceDerivativeMatrix.
*/


/// @brief Convert a local coordinate into a global coordinate.
template <typename Xg, typename Yg, typename X_loc, typename Y_loc>
Coord<Xg, Yg> get_global_coord(Coord<X_loc, Y_loc> const& local_coord)
{
    double const x = ddc::select<X_loc>(local_coord);
    double const y = ddc::select<Y_loc>(local_coord);
    return Coord<Xg, Yg>(x, y);
}

template <typename Xg, typename Yg, typename X_loc, typename Y_loc>
Coord<Xg, Yg> get_global_coord_reverse(
        Coord<X_loc, Y_loc> const& local_coord,
        Coord<X_loc> const x_min,
        Coord<X_loc> const x_max,
        Coord<Y_loc> const y_min,
        Coord<Y_loc> const y_max)
{
    double const x = x_min + x_max - ddc::select<X_loc>(local_coord);
    double const y = y_min + y_max - ddc::select<Y_loc>(local_coord);
    return Coord<Xg, Yg>(x, y);
}

// INITIALISATION OPERATORS ----------------------------------------------------------------------

/// @brief Initialise the function with f(x,y) = cos(2/3*pi*x)sin(y).
template <class Grid1, class Grid2, class Layout>
void initialise_2D_function(DField<IdxRange<Grid1, Grid2>, Kokkos::HostSpace, Layout> function)
{
    ddc::for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        double const xg = ddc::coordinate(Idx<Grid1>(idx));
        double const yg = ddc::coordinate(Idx<Grid2>(idx));
        function(idx) = Kokkos::cos(xg * 2. / 3. * M_PI + 0.25) * Kokkos::sin(yg);
    });
}

template <class Grid1, class Grid2, class Layout>
void initialise_2D_function_reverse(
        DField<IdxRange<Grid1, Grid2>, Kokkos::HostSpace, Layout> function)
{
    IdxRange<Grid1, Grid2> idx_range = get_idx_range(function);
    double const x_max = ddc::coordinate(Idx<Grid1>(idx_range.back()));
    double const y_max = ddc::coordinate(Idx<Grid2>(idx_range.back()));
    double const x_min = ddc::coordinate(Idx<Grid1>(idx_range.front()));
    double const y_min = ddc::coordinate(Idx<Grid2>(idx_range.front()));
    ddc::for_each(idx_range, [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        double const xg = x_min + x_max - ddc::coordinate(Idx<Grid1>(idx));
        double const yg = y_min + y_max - ddc::coordinate(Idx<Grid2>(idx));
        function(idx) = Kokkos::cos(xg * 2. / 3. * M_PI + 0.25) * Kokkos::sin(yg);
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

/// @brief Initialise the y-derivatives of a given DerivField from the global spline.
template <
        class Patch,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type>
void initialise_y_derivatives(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef)
{
    using DerivX = ddc::Deriv<typename Patch::Dim1>;
    using DerivY = ddc::Deriv<typename Patch::Dim2>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    IdxRange<GridX, GridY> idx_range_xy = get_idx_range(function_and_derivs.get_values_field());
    IdxRange<GridX> idx_range_x(idx_range_xy);
    IdxRange<GridY> idx_range_y(idx_range_xy);

    Idx<DerivY, GridY> idx_slice_ymin(Idx<DerivY>(1), idx_range_slice_dy.front());
    Idx<DerivY, GridY> idx_slice_ymax(Idx<DerivY>(1), idx_range_slice_dy.back());

    DField<IdxRange<GridX>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_ymin_extracted
            = function_and_derivs[idx_slice_ymin];
    DField<IdxRange<GridX>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_ymax_extracted
            = function_and_derivs[idx_slice_ymax];

    ddc::for_each(idx_range_x, [&](Idx<GridX> const& idx_par) {
        Idx<GridX, GridY> idx_min(idx_par, idx_range_y.front());
        Idx<GridX, GridY> idx_max(idx_par, idx_range_y.back());
        Coord<Xg, Yg> interface_coord_min(get_global_coord<Xg, Yg>(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(get_global_coord<Xg, Yg>(ddc::coordinate(idx_max)));

        derivs_ymin_extracted(idx_par)
                = evaluator_g.deriv_dim_2(interface_coord_min, const_function_g_coef);
        derivs_ymax_extracted(idx_par)
                = evaluator_g.deriv_dim_2(interface_coord_max, const_function_g_coef);
    });
}

template <
        class Patch,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type>
void initialise_y_derivatives_reversed(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef)
{
    using DerivX = ddc::Deriv<typename Patch::Dim1>;
    using DerivY = ddc::Deriv<typename Patch::Dim2>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    IdxRange<GridX, GridY> idx_range_xy = get_idx_range(function_and_derivs.get_values_field());
    IdxRange<GridX> idx_range_x(idx_range_xy);
    IdxRange<GridY> idx_range_y(idx_range_xy);

    Idx<DerivY, GridY> idx_slice_ymin(Idx<DerivY>(1), idx_range_slice_dy.front());
    Idx<DerivY, GridY> idx_slice_ymax(Idx<DerivY>(1), idx_range_slice_dy.back());

    DField<IdxRange<GridX>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_ymin_extracted
            = function_and_derivs[idx_slice_ymin];
    DField<IdxRange<GridX>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_ymax_extracted
            = function_and_derivs[idx_slice_ymax];

    typename Patch::Coord1 x_min(ddc::coordinate(idx_range_x.front()));
    typename Patch::Coord1 x_max(ddc::coordinate(idx_range_x.back()));
    typename Patch::Coord2 y_min(ddc::coordinate(idx_range_y.front()));
    typename Patch::Coord2 y_max(ddc::coordinate(idx_range_y.back()));
    ddc::for_each(idx_range_x, [&](Idx<GridX> const& idx_par) {
        Idx<GridX, GridY> idx_min(idx_par, idx_range_y.front());
        Idx<GridX, GridY> idx_max(idx_par, idx_range_y.back());
        Coord<Xg, Yg> interface_coord_min(
                get_global_coord_reverse<
                        Xg,
                        Yg>(ddc::coordinate(idx_min), x_min, x_max, y_min, y_max));
        Coord<Xg, Yg> interface_coord_max(
                get_global_coord_reverse<
                        Xg,
                        Yg>(ddc::coordinate(idx_max), x_min, x_max, y_min, y_max));

        derivs_ymin_extracted(idx_par)
                = -evaluator_g.deriv_dim_2(interface_coord_min, const_function_g_coef);
        derivs_ymax_extracted(idx_par)
                = -evaluator_g.deriv_dim_2(interface_coord_max, const_function_g_coef);
    });
}


/// @brief Initialise all the y-derivatives of the given DerivFields from the global spline.
template <
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class... Patches>
void initialise_all_y_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_ranges_slice_dy,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef)
{
    (initialise_y_derivatives<Patches>(
             functions_and_derivs.template get<Patches>(),
             idx_ranges_slice_dy.template get<Patches>(),
             evaluator_g,
             const_function_g_coef),
     ...);
}

/// @brief Initialise the cross-derivatives of a given DerivField from the global spline.
template <
        class Patch,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type>
void initialise_cross_derivatives(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        IdxRange1SliceOnPatch<Patch> const& idx_range_slice_dx,
        IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef,
        double const& xshift)
{
    using DerivX = ddc::Deriv<typename Patch::Dim1>;
    using DerivY = ddc::Deriv<typename Patch::Dim2>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    using IdxdXXdYY = Idx<DerivX, GridX, DerivY, GridY>;

    IdxdXXdYY idx_cross_deriv_min_min(
            Idx<DerivX>(1),
            idx_range_slice_dx.front(),
            Idx<DerivY>(1),
            idx_range_slice_dy.front());
    IdxdXXdYY idx_cross_deriv_max_min(
            Idx<DerivX>(1),
            idx_range_slice_dx.back(),
            Idx<DerivY>(1),
            idx_range_slice_dy.front());
    IdxdXXdYY idx_cross_deriv_min_max(
            Idx<DerivX>(1),
            idx_range_slice_dx.front(),
            Idx<DerivY>(1),
            idx_range_slice_dy.back());
    IdxdXXdYY idx_cross_deriv_max_max(
            Idx<DerivX>(1),
            idx_range_slice_dx.back(),
            Idx<DerivY>(1),
            idx_range_slice_dy.back());

    function_and_derivs(idx_cross_deriv_min_min)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(0 + xshift, 0), const_function_g_coef);
    function_and_derivs(idx_cross_deriv_max_min)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(1 + xshift, 0), const_function_g_coef);
    function_and_derivs(idx_cross_deriv_min_max)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(0 + xshift, 1), const_function_g_coef);
    function_and_derivs(idx_cross_deriv_max_max)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(1 + xshift, 1), const_function_g_coef);
}

template <
        class Patch,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type>
void initialise_cross_derivatives_reverse(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        IdxRange1SliceOnPatch<Patch> const& idx_range_slice_dx,
        IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef,
        double const& xshift)
{
    using DerivX = ddc::Deriv<typename Patch::Dim1>;
    using DerivY = ddc::Deriv<typename Patch::Dim2>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    using IdxdXXdYY = Idx<DerivX, GridX, DerivY, GridY>;

    IdxdXXdYY idx_cross_deriv_min_min(
            Idx<DerivX>(1),
            idx_range_slice_dx.front(),
            Idx<DerivY>(1),
            idx_range_slice_dy.front());
    IdxdXXdYY idx_cross_deriv_max_min(
            Idx<DerivX>(1),
            idx_range_slice_dx.back(),
            Idx<DerivY>(1),
            idx_range_slice_dy.front());
    IdxdXXdYY idx_cross_deriv_min_max(
            Idx<DerivX>(1),
            idx_range_slice_dx.front(),
            Idx<DerivY>(1),
            idx_range_slice_dy.back());
    IdxdXXdYY idx_cross_deriv_max_max(
            Idx<DerivX>(1),
            idx_range_slice_dx.back(),
            Idx<DerivY>(1),
            idx_range_slice_dy.back());

    function_and_derivs(idx_cross_deriv_max_max)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(0 + xshift, 0), const_function_g_coef);
    function_and_derivs(idx_cross_deriv_min_max)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(1 + xshift, 0), const_function_g_coef);
    function_and_derivs(idx_cross_deriv_max_min)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(0 + xshift, 1), const_function_g_coef);
    function_and_derivs(idx_cross_deriv_min_min)
            = evaluator_g.deriv_1_and_2(Coord<Xg, Yg>(1 + xshift, 1), const_function_g_coef);
}

/// @brief Initialise all the cross-derivatives of the given DerivFields from the global spline.
template <
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Patch1,
        class Patch2,
        class Patch3>
void initialise_all_cross_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patch1, Patch2, Patch3>& functions_and_derivs,
        MultipatchType<IdxRange1SliceOnPatch, Patch1, Patch2, Patch3> const& idx_ranges_slice_dx,
        MultipatchType<IdxRange2SliceOnPatch, Patch1, Patch2, Patch3> const& idx_ranges_slice_dy,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef)
{
    initialise_cross_derivatives<Patch1>(
            functions_and_derivs.template get<Patch1>(),
            idx_ranges_slice_dx.template get<Patch1>(),
            idx_ranges_slice_dy.template get<Patch1>(),
            evaluator_g,
            const_function_g_coef,
            0);
    initialise_cross_derivatives<Patch2>(
            functions_and_derivs.template get<Patch2>(),
            idx_ranges_slice_dx.template get<Patch2>(),
            idx_ranges_slice_dy.template get<Patch2>(),
            evaluator_g,
            const_function_g_coef,
            1);
    initialise_cross_derivatives<Patch3>(
            functions_and_derivs.template get<Patch3>(),
            idx_ranges_slice_dx.template get<Patch3>(),
            idx_ranges_slice_dy.template get<Patch3>(),
            evaluator_g,
            const_function_g_coef,
            2);
}

// CHECK OPERATORS -------------------------------------------------------------------------------

/** @brief Check that the local grids and the equivalent global grid match together for a 
     * given patch. The integers x_shift and y_shift correspond to a shift in the indices to match
     * with the correct index on the global grid. 
     */
template <class Patch, class GridXg, class GridYg>
void check_interpolation_grids(
        typename Patch::IdxRange12 const& idx_range,
        int const x_shift,
        int const y_shift = 0)
{
    using Grid1 = typename Patch::Grid1;
    using Grid2 = typename Patch::Grid2;
    ddc::for_each(idx_range, [&](typename Patch::Idx12 const& idx) {
        IdxStep<Grid1> idx_x(Idx<Grid1>(idx) - IdxRange<Grid1>(idx_range).front());
        IdxStep<Grid2> idx_y(Idx<Grid2>(idx) - IdxRange<Grid2>(idx_range).front());
        Idx<GridXg, GridYg> idx_g(idx_x.value() + x_shift, idx_y.value() + y_shift);
        EXPECT_NEAR(ddc::coordinate(Idx<Grid1>(idx)), ddc::coordinate(Idx<GridXg>(idx_g)), 1e-15);
        EXPECT_NEAR(ddc::coordinate(Idx<Grid2>(idx)), ddc::coordinate(Idx<GridYg>(idx_g)), 1e-15);
    });
}

template <class Patch, class GridXg, class GridYg>
void check_interpolation_grids_reverse(
        typename Patch::IdxRange12 const& idx_range,
        int const x_shift,
        int const y_shift = 0)
{
    using Grid1 = typename Patch::Grid1;
    using Grid2 = typename Patch::Grid2;
    using Dim1 = typename Patch::Dim1;
    using Dim2 = typename Patch::Dim2;
    Coord<Dim1> coord_x_max = ddc::coordinate(Idx<Grid1>(idx_range.back()));
    Coord<Dim1> coord_x_min = ddc::coordinate(Idx<Grid1>(idx_range.front()));
    Coord<Dim2> coord_y_max = ddc::coordinate(Idx<Grid2>(idx_range.back()));
    Coord<Dim2> coord_y_min = ddc::coordinate(Idx<Grid2>(idx_range.front()));
    ddc::for_each(idx_range, [&](typename Patch::Idx12 const& idx) {
        IdxStep<Grid1> idx_x(IdxRange<Grid1>(idx_range).back() - Idx<Grid1>(idx));
        IdxStep<Grid2> idx_y(IdxRange<Grid2>(idx_range).back() - Idx<Grid2>(idx));
        Idx<GridXg, GridYg> idx_g(idx_x.value() + x_shift, idx_y.value() + y_shift);
        EXPECT_NEAR(
                coord_x_max + coord_x_min - ddc::coordinate(Idx<Grid1>(idx)),
                ddc::coordinate(Idx<GridXg>(idx_g)),
                1e-15);
        EXPECT_NEAR(
                coord_y_max + coord_y_min - ddc::coordinate(Idx<Grid2>(idx)),
                ddc::coordinate(Idx<GridYg>(idx_g)),
                1e-15);
    });
}


/** @brief Check agreement between the computed x-derivatives and the global x-derivatives at 
     * the interfaces for a given patch.
     */
template <
        class Patch,
        class ReversedPatchSeq,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2>
void check_x_derivatives(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        typename Patch::IdxRange1 const& idx_range_perp,
        IdxRange1SliceOnPatch<Patch> const& idx_range_slice_dx)
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    constexpr bool is_reversed_patch = ddc::in_tags_v<Patch, ReversedPatchSeq>;
    constexpr int reversing_sign = !is_reversed_patch - is_reversed_patch;

    using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
    Idx<DerivX> idx_deriv(1);

    Idx<DerivX, typename Patch::Grid1> idx_slice_min(idx_deriv, idx_range_slice_dx.front());
    Idx<DerivX, typename Patch::Grid1> idx_slice_max(idx_deriv, idx_range_slice_dx.back());

    DField<IdxRange<typename Patch::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_xmin_extracted = function_and_derivs[idx_slice_min];
    DField<IdxRange<typename Patch::Grid2>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_xmax_extracted = function_and_derivs[idx_slice_max];

    typename Patch::IdxRange2 idx_range_par(get_idx_range(derivs_xmin_extracted));

    typename Patch::Coord1 x_min(ddc::coordinate(idx_range_perp.front()));
    typename Patch::Coord1 x_max(ddc::coordinate(idx_range_perp.back()));
    typename Patch::Coord2 y_min(ddc::coordinate(idx_range_par.front()));
    typename Patch::Coord2 y_max(ddc::coordinate(idx_range_par.back()));
    ddc::for_each(idx_range_par, [&](typename Patch::Idx2 const& idx_par) {
        typename Patch::Idx12 idx_min(idx_range_perp.front(), idx_par);
        typename Patch::Idx12 idx_max(idx_range_perp.back(), idx_par);
        Coord<Xg, Yg> interface_coord_min;
        Coord<Xg, Yg> interface_coord_max;
        if (is_reversed_patch) {
            interface_coord_min = get_global_coord_reverse<
                    Xg,
                    Yg>(ddc::coordinate(idx_min), x_min, x_max, y_min, y_max);
            interface_coord_max = get_global_coord_reverse<
                    Xg,
                    Yg>(ddc::coordinate(idx_max), x_min, x_max, y_min, y_max);
        } else {
            interface_coord_min = get_global_coord<Xg, Yg>(ddc::coordinate(idx_min));
            interface_coord_max = get_global_coord<Xg, Yg>(ddc::coordinate(idx_max));
        }

        double const global_deriv_min
                = reversing_sign * evaluator_g.deriv_dim_1(interface_coord_min, function_g_coef);
        double const global_deriv_max
                = reversing_sign * evaluator_g.deriv_dim_1(interface_coord_max, function_g_coef);

        EXPECT_NEAR(derivs_xmin_extracted(idx_par), global_deriv_min, 5e-14);
        EXPECT_NEAR(derivs_xmax_extracted(idx_par), global_deriv_max, 5e-14);
    });
}

/** @brief Check agreement between the computed x-derivatives and the global x-derivatives at 
     * the interfaces. 
     */
template <
        class ReversedPatchSeq,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class... Patches>
void check_all_x_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_range_slices_dx)
{
    (check_x_derivatives<Patches, ReversedPatchSeq>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             typename Patches::IdxRange1(idx_ranges.template get<Patches>()),
             idx_range_slices_dx.template get<Patches>()),
     ...);
}



/** @brief Check agreement between the computed y-derivatives and the global y-derivatives at 
     * the interfaces for a given patch.
     */
template <
        class Patch,
        class ReversedPatchSeq,
        class PatchSeqMin,
        class PatchSeqMax,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2>
void check_y_derivatives(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        typename Patch::IdxRange2 const& idx_range_perp,
        IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy)
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    constexpr bool is_reversed_patch = ddc::in_tags_v<Patch, ReversedPatchSeq>;
    constexpr int reversing_sign = !is_reversed_patch - is_reversed_patch;

    using DerivY = typename ddc::Deriv<typename Patch::Dim2>;
    Idx<DerivY> idx_deriv(1);

    Idx<DerivY, typename Patch::Grid2> idx_slice_min(idx_deriv, idx_range_slice_dy.front());
    Idx<DerivY, typename Patch::Grid2> idx_slice_max(idx_deriv, idx_range_slice_dy.back());

    DField<IdxRange<typename Patch::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymin_extracted = function_and_derivs[idx_slice_min];
    DField<IdxRange<typename Patch::Grid1>, Kokkos::HostSpace, Kokkos::layout_stride>
            derivs_ymax_extracted = function_and_derivs[idx_slice_max];

    typename Patch::IdxRange1 idx_range_par(get_idx_range(derivs_ymin_extracted));
    typename Patch::Idx2 idx_y_min = idx_range_perp.front();
    typename Patch::Idx2 idx_y_max = idx_range_perp.back();

    typename Patch::Coord1 x_min(ddc::coordinate(idx_range_par.front()));
    typename Patch::Coord1 x_max(ddc::coordinate(idx_range_par.back()));
    typename Patch::Coord2 y_min(ddc::coordinate(idx_range_perp.front()));
    typename Patch::Coord2 y_max(ddc::coordinate(idx_range_perp.back()));
    ddc::for_each(idx_range_par, [&](typename Patch::Idx1 const& idx) {
        typename Patch::Idx12 idx_min(idx, idx_y_min);
        typename Patch::Idx12 idx_max(idx, idx_y_max);
        Coord<Xg, Yg> interface_coord_min;
        Coord<Xg, Yg> interface_coord_max;

        if (is_reversed_patch) {
            interface_coord_min = get_global_coord_reverse<
                    Xg,
                    Yg>(ddc::coordinate(idx_min), x_min, x_max, y_min, y_max);
            interface_coord_max = get_global_coord_reverse<
                    Xg,
                    Yg>(ddc::coordinate(idx_max), x_min, x_max, y_min, y_max);
        } else {
            interface_coord_min = get_global_coord<Xg, Yg>(ddc::coordinate(idx_min));
            interface_coord_max = get_global_coord<Xg, Yg>(ddc::coordinate(idx_max));
        }

        double const global_deriv_min
                = reversing_sign * evaluator_g.deriv_dim_2(interface_coord_min, function_g_coef);
        double const global_deriv_max
                = reversing_sign * evaluator_g.deriv_dim_2(interface_coord_max, function_g_coef);

        // For Patches in PatchSeqMin, we defined ddc::BoundCond::GREVILLE the local lower Y-boundary,
        // we don't need the y-derivatives for y = ymin. Their values are not checked.
        if constexpr (!ddc::in_tags_v<Patch, PatchSeqMin>) {
            EXPECT_NEAR(derivs_ymin_extracted(idx), global_deriv_min, 5e-14);
        }
        // For Patches in PatchSeqMax, we defined ddc::BoundCond::GREVILLE the local lower Y-boundary,
        // we don't need the y-derivatives for y = ymax. Their values are not checked.
        if constexpr (!ddc::in_tags_v<Patch, PatchSeqMax>) {
            EXPECT_NEAR(derivs_ymax_extracted(idx), global_deriv_max, 5e-14);
        }
    });
}

/** @brief Check agreement between the computed y-derivatives and the global y-derivatives at 
     * the interfaces.
     */
template <
        class ReversedPatchSeq,
        class PatchSeqMin,
        class PatchSeqMax,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class... Patches>
void check_all_y_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_range_slices_dy)
{
    (check_y_derivatives<Patches, ReversedPatchSeq, PatchSeqMin, PatchSeqMax>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             typename Patches::IdxRange2(idx_ranges.template get<Patches>()),
             idx_range_slices_dy.template get<Patches>()),
     ...);
}


/** @brief Check agreement between the computed cross-derivatives and the global cross-derivatives at 
     * the interfaces for a given patch.
     */
template <
        class Patch,
        class ReversedPatchSeq,
        class PatchSeqMin,
        class PatchSeqMax,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2>
void check_xy_derivatives(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        typename Patch::IdxRange12 const& idx_range,
        IdxRange1SliceOnPatch<Patch> const& idx_range_slice_dx,
        IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy)
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    constexpr bool is_reversed_patch = ddc::in_tags_v<Patch, ReversedPatchSeq>;
    constexpr int reversing_sign = !is_reversed_patch - is_reversed_patch;

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

    Coord<Xg, Yg> interface_coord_min_min(get_global_coord<Xg, Yg>(ddc::coordinate(idx_min_min)));
    Coord<Xg, Yg> interface_coord_max_min(get_global_coord<Xg, Yg>(ddc::coordinate(idx_max_min)));
    Coord<Xg, Yg> interface_coord_min_max(get_global_coord<Xg, Yg>(ddc::coordinate(idx_min_max)));
    Coord<Xg, Yg> interface_coord_max_max(get_global_coord<Xg, Yg>(ddc::coordinate(idx_max_max)));

    double global_deriv_min_min;
    double global_deriv_max_min;
    double global_deriv_min_max;
    double global_deriv_max_max;
    if (is_reversed_patch) {
        global_deriv_min_min = evaluator_g.deriv_1_and_2(interface_coord_max_max, function_g_coef);
        global_deriv_max_min = evaluator_g.deriv_1_and_2(interface_coord_min_max, function_g_coef);
        global_deriv_min_max = evaluator_g.deriv_1_and_2(interface_coord_max_min, function_g_coef);
        global_deriv_max_max = evaluator_g.deriv_1_and_2(interface_coord_min_min, function_g_coef);
    } else {
        global_deriv_min_min = evaluator_g.deriv_1_and_2(interface_coord_min_min, function_g_coef);
        global_deriv_max_min = evaluator_g.deriv_1_and_2(interface_coord_max_min, function_g_coef);
        global_deriv_min_max = evaluator_g.deriv_1_and_2(interface_coord_min_max, function_g_coef);
        global_deriv_max_max = evaluator_g.deriv_1_and_2(interface_coord_max_max, function_g_coef);
    }

    // For Patches in PatchSeqMin, we defined ddc::BoundCond::GREVILLE the local lower Y-boundary,
    // we don't need the cross-derivatives for y = ymin. Their value is not checked.
    if constexpr (!ddc::in_tags_v<Patch, PatchSeqMin>) {
        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_min_min), global_deriv_min_min, 5e-13);
        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_max_min), global_deriv_max_min, 5e-13);
    }
    // For Patches in PatchSeqMax, we defined ddc::BoundCond::GREVILLE the local upper Y-boundary,
    // we don't need the cross-derivatives for y = ymax. Their value is not checked.
    if constexpr (!ddc::in_tags_v<Patch, PatchSeqMax>) {
        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_min_max), global_deriv_min_max, 5e-13);
        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_max_max), global_deriv_max_max, 5e-13);
    }
}

/** @brief Check agreement between the computed cross-derivatives and the global cross-derivatives at 
     * the interfaces.
     */
template <
        class ReversedPatchSeq,
        class PatchSeqMin,
        class PatchSeqMax,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class... Patches>
void check_all_xy_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_range_slices_dx,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_range_slices_dy)
{
    (check_xy_derivatives<Patches, ReversedPatchSeq, PatchSeqMin, PatchSeqMax>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             idx_ranges.template get<Patches>(),
             idx_range_slices_dx.template get<Patches>(),
             idx_range_slices_dy.template get<Patches>()),
     ...);
}


/** @brief Check agreement between the local splines defined with the computed derivatives 
     * and the global spline for a given patch.  
     */
template <
        class Patch,
        class ReversedPatchSeq,
        class PatchSeqMin,
        class PatchSeqMax,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2>
void check_spline_representation_agreement(
        typename Patch::IdxRange12 const& idx_range_xy,
        IdxRange1SliceOnPatch<Patch> const& idx_range_slice_x,
        IdxRange2SliceOnPatch<Patch> const& idx_range_slice_y,
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
{
    using HostExecSpace = Kokkos::DefaultHostExecutionSpace;

    using DimX = typename Patch::Dim1;
    using DimY = typename Patch::Dim2;
    using DerivX = typename ddc::Deriv<DimX>;
    using DerivY = typename ddc::Deriv<DimY>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    constexpr bool is_reversed_patch = ddc::in_tags_v<Patch, ReversedPatchSeq>;

    const ddc::BoundCond BoundCondXmin = ddc::BoundCond::HERMITE;
    const ddc::BoundCond BoundCondXmax = ddc::BoundCond::HERMITE;

    // For Patches in PatchSeqMin, the local lower Y-boundary is with the outside.
    const ddc::BoundCond BoundCondYmin = ddc::in_tags_v<Patch, PatchSeqMin>
                                                 ? ddc::BoundCond::GREVILLE
                                                 : ddc::BoundCond::HERMITE;
    // For Patches in PatchSeqMax, the local upper Y-boundary is with the outside.
    const ddc::BoundCond BoundCondYmax = ddc::in_tags_v<Patch, PatchSeqMax>
                                                 ? ddc::BoundCond::GREVILLE
                                                 : ddc::BoundCond::HERMITE;

    // Build the local spline representation -------------------------------------------------
    ddc::SplineBuilder2D<
            HostExecSpace,
            typename HostExecSpace::memory_space,
            typename Patch::BSplines1,
            typename Patch::BSplines2,
            typename Patch::Grid1,
            typename Patch::Grid2,
            BoundCondXmin,
            BoundCondXmax,
            BoundCondYmin,
            BoundCondYmax,
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

    DField<IdxRange<DerivX, GridY>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_xmin_extracted
            = function_and_derivs[idx_range_deriv_xmin][idx_slice_xmin];
    DField<IdxRange<DerivX, GridY>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_xmax_extracted
            = function_and_derivs[idx_range_deriv_xmax][idx_slice_xmax];

    DField<IdxRange<GridX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_ymin_extracted
            = function_and_derivs[idx_range_deriv_ymin][idx_slice_ymin];
    DField<IdxRange<GridX, DerivY>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_ymax_extracted
            = function_and_derivs[idx_range_deriv_ymax][idx_slice_ymax];

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

    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_min_min_alloc(idx_range_deriv_x_deriv_y);
    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_max_min_alloc(idx_range_deriv_x_deriv_y);
    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_min_max_alloc(idx_range_deriv_x_deriv_y);
    host_t<DFieldMem<IdxRange<DerivX, DerivY>>> deriv_xy_max_max_alloc(idx_range_deriv_x_deriv_y);

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


    // If the boundary is not a ddc::BoundCond::HERMITE, we don't use derivatives.
    if constexpr (
            (BoundCondYmin == ddc::BoundCond::HERMITE)
            && (BoundCondYmax == ddc::BoundCond::HERMITE)) {
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
    } else if constexpr (
            (BoundCondYmin == ddc::BoundCond::HERMITE)
            && !(BoundCondYmax == ddc::BoundCond::HERMITE)) {
        builder(function_coef,
                get_const_field(function),
                std::optional(get_const_field(deriv_xmin)),
                std::optional(get_const_field(deriv_xmax)),
                std::optional(get_const_field(deriv_ymin)),
                std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
                std::optional(get_const_field(deriv_xy_min_min)),
                std::optional(get_const_field(deriv_xy_max_min)),
                std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt});
    } else if constexpr (
            !(BoundCondYmin == ddc::BoundCond::HERMITE)
            && (BoundCondYmax == ddc::BoundCond::HERMITE)) {
        builder(function_coef,
                get_const_field(function),
                std::optional(get_const_field(deriv_xmin)),
                std::optional(get_const_field(deriv_xmax)),
                std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
                std::optional(get_const_field(deriv_ymax)),
                std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                std::optional(get_const_field(deriv_xy_min_max)),
                std::optional(get_const_field(deriv_xy_max_max)));
    } else {
        builder(function_coef,
                get_const_field(function),
                std::optional(get_const_field(deriv_xmin)),
                std::optional(get_const_field(deriv_xmax)),
                std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
                std::optional<ConstDeriv2_OnPatch_2D_host<Patch>> {std::nullopt},
                std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt},
                std::optional<ConstDeriv12_OnPatch_2D_host<Patch>> {std::nullopt});
    }


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
        Coord<Xg, Yg> eval_point_g;
        if (is_reversed_patch) {
            eval_point_g = get_global_coord_reverse<
                    Xg,
                    Yg>(eval_points(idx), x_min, x_max, y_min, y_max);
        } else {
            eval_point_g = (get_global_coord<Xg, Yg>(eval_points(idx)));
        }

        double local_spline = evaluator(eval_points(idx), get_const_field(function_coef));
        double global_spline = evaluator_g(eval_point_g, get_const_field(function_g_coef));

        EXPECT_NEAR(local_spline, global_spline, 1e-14);
    });
}

/** @brief Check agreement between the local splines defined with the computed derivatives 
     * and the global spline. 
     */
template <
        class ReversedPatchSeq,
        class PatchSeqMin,
        class PatchSeqMax,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class... Patches>
void check_all_spline_representation_agreement(
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_range_slices_x,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_range_slices_y,
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef)
{
    (check_spline_representation_agreement<Patches, ReversedPatchSeq, PatchSeqMin, PatchSeqMax>(
             idx_ranges.template get<Patches>(),
             idx_range_slices_x.template get<Patches>(),
             idx_range_slices_y.template get<Patches>(),
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             get_const_field(function_g_coef)),
     ...);
};