
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "multipatch_field.hpp"
#include "spline_builder_deriv_field.hpp"
#include "types.hpp"

/*
    This file defined useful operators for the test on the InterfaceExactDerivativeMatrix.
*/

/// @brief Coordinate transformation from the patches to the global/physical domain.
template <typename Xg, typename Yg, typename X_loc, typename Y_loc>
struct CoordTransform
{
    using Dim1Local = X_loc;
    using Dim2Local = Y_loc;
    using Dim1Global = Xg;
    using Dim2Global = Yg;

    const bool m_is_reverse_x;
    const bool m_is_reverse_y;
    const bool m_are_exchange_x_y;

    const Coord<X_loc> m_x_min;
    const Coord<X_loc> m_x_max;
    const Coord<Y_loc> m_y_min;
    const Coord<Y_loc> m_y_max;

    /**
     * @brief Instantiate the coordinate transformator.
     * By default, the transformation is identity.
     */
    CoordTransform(
            bool is_reverse_x = false,
            bool is_reverse_y = false,
            bool are_exchange_x_y = false,
            Coord<X_loc> x_min = Coord<X_loc>(0),
            Coord<X_loc> x_max = Coord<X_loc>(1),
            Coord<Y_loc> y_min = Coord<Y_loc>(0),
            Coord<Y_loc> y_max = Coord<Y_loc>(1))
        : m_is_reverse_x(is_reverse_x)
        , m_is_reverse_y(is_reverse_y)
        , m_are_exchange_x_y(are_exchange_x_y)
        , m_x_min(x_min)
        , m_x_max(x_max)
        , m_y_min(y_min)
        , m_y_max(y_max)
    {
    }

    Coord<Xg, Yg> get_global_coord(Coord<X_loc, Y_loc> const& local_coord) const
    {
        double x_loc = ddc::select<X_loc>(local_coord);
        double y_loc = ddc::select<Y_loc>(local_coord);

        if (m_is_reverse_x) {
            x_loc = m_x_min + m_x_max - x_loc;
        }
        if (m_is_reverse_y) {
            y_loc = m_y_min + m_y_max - y_loc;
        }
        if (m_are_exchange_x_y) {
            return Coord<Xg, Yg>(y_loc, x_loc);
        } else {
            return Coord<Xg, Yg>(x_loc, y_loc);
        }
    }
};

// -----------------------------------------------------------------------------------------------
// INITIALISATION OPERATORS ----------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------

/// @brief Initialise the function with f(x,y) = cos(2/3*pi*x)sin(y).
template <class Grid1, class Grid2, class CoordTransformType, class Layout>
void initialise_2D_function(
        DField<IdxRange<Grid1, Grid2>, Kokkos::HostSpace, Layout> function,
        CoordTransformType const& coord_transform = CoordTransformType())
{
    using Dim1 = typename CoordTransformType::Dim1Local;
    using Dim2 = typename CoordTransformType::Dim2Local;
    using Xg = typename CoordTransformType::Dim1Global;
    using Yg = typename CoordTransformType::Dim2Global;
    ddc::for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        Coord<Dim1, Dim2> local_coord = ddc::coordinate(idx);
        Coord<Xg, Yg> equiv_global_coord = coord_transform.get_global_coord(local_coord);
        double const xg = ddc::get<Xg>(equiv_global_coord);
        double const yg = ddc::get<Yg>(equiv_global_coord);
        function(idx) = Kokkos::cos(xg * 2. / 3. * M_PI + 0.25) * Kokkos::sin(yg);
    });
}

/// @brief Initialise all the functions defined on the patches.
template <class Xg, class Yg, class... Patches>
void initialise_all_functions(
        MultipatchField<DerivFieldOnPatch_host, Patches...> const& functions_and_derivs,
        std::tuple<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>...>
                coord_transforms)
{
    (initialise_2D_function<typename Patches::Grid1, typename Patches::Grid2>(
             (functions_and_derivs.template get<Patches>()).get_values_field(),
             std::get<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>>(
                     coord_transforms)),
     ...);
}


/// @brief Initialise the y-derivatives of a given DerivField from the global spline.
template <
        class Patch,
        class CoordTransformType,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type>
void initialise_y_derivatives(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef,
        CoordTransformType const& coord_transform = CoordTransformType())
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
        Coord<Xg, Yg> interface_coord_min(
                coord_transform.get_global_coord(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(
                coord_transform.get_global_coord(ddc::coordinate(idx_max)));

        int const sign = !coord_transform.m_is_reverse_y - coord_transform.m_is_reverse_y;

        if (coord_transform.m_are_exchange_x_y) {
            derivs_ymin_extracted(idx_par)
                    = sign * evaluator_g.deriv_dim_1(interface_coord_min, const_function_g_coef);
            derivs_ymax_extracted(idx_par)
                    = sign * evaluator_g.deriv_dim_1(interface_coord_max, const_function_g_coef);
        } else {
            derivs_ymin_extracted(idx_par)
                    = sign * evaluator_g.deriv_dim_2(interface_coord_min, const_function_g_coef);
            derivs_ymax_extracted(idx_par)
                    = sign * evaluator_g.deriv_dim_2(interface_coord_max, const_function_g_coef);
        }
    });
}


/// @brief Initialise the x-derivatives of a given DerivField from the global spline.
template <
        class Patch,
        class CoordTransformType,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type>
void initialise_x_derivatives(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        IdxRange1SliceOnPatch<Patch> const& idx_range_slice_dx,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef,
        CoordTransformType const& coord_transform = CoordTransformType())
{
    using DerivX = ddc::Deriv<typename Patch::Dim1>;
    using DerivY = ddc::Deriv<typename Patch::Dim2>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    IdxRange<GridX, GridY> idx_range_xy = get_idx_range(function_and_derivs.get_values_field());
    IdxRange<GridX> idx_range_x(idx_range_xy);
    IdxRange<GridY> idx_range_y(idx_range_xy);

    Idx<DerivX, GridX> idx_slice_xmin(Idx<DerivX>(1), idx_range_slice_dx.front());
    Idx<DerivX, GridX> idx_slice_xmax(Idx<DerivX>(1), idx_range_slice_dx.back());

    DField<IdxRange<GridY>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_xmin_extracted
            = function_and_derivs[idx_slice_xmin];
    DField<IdxRange<GridY>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_xmax_extracted
            = function_and_derivs[idx_slice_xmax];

    ddc::for_each(idx_range_y, [&](Idx<GridY> const& idx_par) {
        Idx<GridX, GridY> idx_min(idx_range_x.front(), idx_par);
        Idx<GridX, GridY> idx_max(idx_range_x.back(), idx_par);
        Coord<Xg, Yg> interface_coord_min(
                coord_transform.get_global_coord(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(
                coord_transform.get_global_coord(ddc::coordinate(idx_max)));

        int const sign = !coord_transform.m_is_reverse_x - coord_transform.m_is_reverse_x;

        if (coord_transform.m_are_exchange_x_y) {
            derivs_xmin_extracted(idx_par)
                    = sign * evaluator_g.deriv_dim_2(interface_coord_min, const_function_g_coef);
            derivs_xmax_extracted(idx_par)
                    = sign * evaluator_g.deriv_dim_2(interface_coord_max, const_function_g_coef);
        } else {
            derivs_xmin_extracted(idx_par)
                    = sign * evaluator_g.deriv_dim_1(interface_coord_min, const_function_g_coef);
            derivs_xmax_extracted(idx_par)
                    = sign * evaluator_g.deriv_dim_1(interface_coord_max, const_function_g_coef);
        }
    });
}


/// @brief Initialise the cross-derivatives of a given DerivField from the global spline.
template <
        class Patch,
        class CoordTransformType,
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
        double const& xshift,
        CoordTransformType const& coord_transform = CoordTransformType())
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

    Coord<Xg, Yg> coord_min_min = coord_transform.get_global_coord(
            typename Patch::Coord12(coord_transform.m_x_min, coord_transform.m_y_min));
    Coord<Xg, Yg> coord_max_min = coord_transform.get_global_coord(
            typename Patch::Coord12(coord_transform.m_x_max, coord_transform.m_y_min));
    Coord<Xg, Yg> coord_min_max = coord_transform.get_global_coord(
            typename Patch::Coord12(coord_transform.m_x_min, coord_transform.m_y_max));
    Coord<Xg, Yg> coord_max_max = coord_transform.get_global_coord(
            typename Patch::Coord12(coord_transform.m_x_max, coord_transform.m_y_max));

    int const sign_x = !coord_transform.m_is_reverse_x - coord_transform.m_is_reverse_x;
    int const sign_y = !coord_transform.m_is_reverse_y - coord_transform.m_is_reverse_y;
    int const sign_xy = sign_x * sign_y;

    function_and_derivs(idx_cross_deriv_min_min)
            = sign_xy * evaluator_g.deriv_1_and_2(coord_min_min, const_function_g_coef);
    function_and_derivs(idx_cross_deriv_max_min)
            = sign_xy * evaluator_g.deriv_1_and_2(coord_max_min, const_function_g_coef);
    function_and_derivs(idx_cross_deriv_min_max)
            = sign_xy * evaluator_g.deriv_1_and_2(coord_min_max, const_function_g_coef);
    function_and_derivs(idx_cross_deriv_max_max)
            = sign_xy * evaluator_g.deriv_1_and_2(coord_max_max, const_function_g_coef);
}

/// @brief Initialise all the cross-derivatives of the given DerivFields from the global spline.
template <
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type,
        class... Patches>
void initialise_all_cross_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_ranges_slice_dx,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_ranges_slice_dy,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef,
        std::tuple<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>...>
                coord_transforms)
{
    using PatchSeq = ddc::detail::TypeSeq<Patches...>;
    (initialise_cross_derivatives<Patches>(
             functions_and_derivs.template get<Patches>(),
             idx_ranges_slice_dx.template get<Patches>(),
             idx_ranges_slice_dy.template get<Patches>(),
             evaluator_g,
             const_function_g_coef,
             ddc::type_seq_rank_v<Patches, PatchSeq>,
             std::get<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>>(
                     coord_transforms)),
     ...);
}

// -----------------------------------------------------------------------------------------------
// CHECK OPERATORS -------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------

/** 
 * @brief Check that the local grids and the equivalent global grid match together for a 
 * given patch. The integers x_shift and y_shift correspond to a shift in the indices to match
 * with the correct index on the global grid. 
 */
template <
        class Patch,
        class GridXg,
        class GridYg,
        class CoordTransformType = CoordTransform<
                typename GridXg::continuous_dimension_type,
                typename GridYg::continuous_dimension_type,
                typename Patch::Dim1,
                typename Patch::Dim2>>
void check_interpolation_grids(
        typename Patch::IdxRange12 const& idx_range,
        int const x_shift,
        int const y_shift,
        CoordTransformType const coord_transform = CoordTransformType())
{
    using Grid1 = typename Patch::Grid1;
    using Grid2 = typename Patch::Grid2;
    using Xg = typename GridXg::continuous_dimension_type;
    using Yg = typename GridYg::continuous_dimension_type;
    using Dim1 = typename Patch::Dim1;
    using Dim2 = typename Patch::Dim2;
    ddc::for_each(idx_range, [&](typename Patch::Idx12 const& idx) {
        IdxStep<Grid1> idx_x = coord_transform.m_is_reverse_x
                                       ? IdxRange<Grid1>(idx_range).back() - Idx<Grid1>(idx)
                                       : Idx<Grid1>(idx) - IdxRange<Grid1>(idx_range).front();
        IdxStep<Grid2> idx_y = coord_transform.m_is_reverse_y
                                       ? IdxRange<Grid2>(idx_range).back() - Idx<Grid2>(idx)
                                       : Idx<Grid2>(idx) - IdxRange<Grid2>(idx_range).front();
        Idx<GridXg, GridYg> idx_g
                = coord_transform.m_are_exchange_x_y
                          ? Idx<GridXg, GridYg>(idx_y.value() + x_shift, idx_x.value() + y_shift)
                          : Idx<GridXg, GridYg>(idx_x.value() + x_shift, idx_y.value() + y_shift);

        Coord<Dim1, Dim2> local_coord = ddc::coordinate(idx);
        Coord<Xg, Yg> equiv_global_coord = coord_transform.get_global_coord(local_coord);
        Coord<Xg, Yg> global_coord = ddc::coordinate(idx_g);
        EXPECT_NEAR(ddc::get<Xg>(equiv_global_coord), ddc::get<Xg>(global_coord), 1e-15);
        EXPECT_NEAR(ddc::get<Yg>(equiv_global_coord), ddc::get<Yg>(global_coord), 1e-15);
    });
}


/** 
 * @brief Check agreement between the computed x-derivatives and the global x-derivatives at 
 * the interfaces for a given patch.
 */
template <
        class Patch,
        class CoordTransformType,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2>
void check_x_derivatives(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        typename Patch::IdxRange1 const& idx_range_perp,
        IdxRange1SliceOnPatch<Patch> const& idx_range_slice_dx,
        CoordTransformType const& coord_transform = CoordTransformType())
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    const int sign = !coord_transform.m_is_reverse_x - coord_transform.m_is_reverse_x;

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
        Coord<Xg, Yg> interface_coord_min(
                coord_transform.get_global_coord(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(
                coord_transform.get_global_coord(ddc::coordinate(idx_max)));

        double global_deriv_min;
        double global_deriv_max;
        if (coord_transform.m_are_exchange_x_y) {
            global_deriv_min = sign * evaluator_g.deriv_dim_2(interface_coord_min, function_g_coef);
            global_deriv_max = sign * evaluator_g.deriv_dim_2(interface_coord_max, function_g_coef);

        } else {
            global_deriv_min = sign * evaluator_g.deriv_dim_1(interface_coord_min, function_g_coef);
            global_deriv_max = sign * evaluator_g.deriv_dim_1(interface_coord_max, function_g_coef);
        }
        EXPECT_NEAR(derivs_xmin_extracted(idx_par), global_deriv_min, 5e-14);
        EXPECT_NEAR(derivs_xmax_extracted(idx_par), global_deriv_max, 5e-14);
    });
}

/** 
 * @brief Check agreement between the computed x-derivatives and the global x-derivatives at 
 * the interfaces. 
 */
template <
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type,
        class... Patches>
void check_all_x_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_range_slices_dx,
        std::tuple<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>...>
                coord_transforms)
{
    (check_x_derivatives<Patches>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             typename Patches::IdxRange1(idx_ranges.template get<Patches>()),
             idx_range_slices_dx.template get<Patches>(),
             std::get<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>>(
                     coord_transforms)),
     ...);
}


/** 
 * @brief Check agreement between the computed y-derivatives and the global y-derivatives at 
 * the interfaces for a given patch.
 */
template <
        class Patch,
        class PatchSeqMin,
        class PatchSeqMax,
        class CoordTransformType,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2>
void check_y_derivatives(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        typename Patch::IdxRange2 const& idx_range_perp,
        IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy,
        CoordTransformType const& coord_transform = CoordTransformType())
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    const int sign = !coord_transform.m_is_reverse_y - coord_transform.m_is_reverse_y;

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

    ddc::for_each(idx_range_par, [&](typename Patch::Idx1 const& idx) {
        typename Patch::Idx12 idx_min(idx, idx_y_min);
        typename Patch::Idx12 idx_max(idx, idx_y_max);
        Coord<Xg, Yg> interface_coord_min(
                coord_transform.get_global_coord(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(
                coord_transform.get_global_coord(ddc::coordinate(idx_max)));

        double global_deriv_min;
        double global_deriv_max;
        if (coord_transform.m_are_exchange_x_y) {
            global_deriv_min = sign * evaluator_g.deriv_dim_1(interface_coord_min, function_g_coef);
            global_deriv_max = sign * evaluator_g.deriv_dim_1(interface_coord_max, function_g_coef);

        } else {
            global_deriv_min = sign * evaluator_g.deriv_dim_2(interface_coord_min, function_g_coef);
            global_deriv_max = sign * evaluator_g.deriv_dim_2(interface_coord_max, function_g_coef);
        }

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

/** 
 * @brief Check agreement between the computed y-derivatives and the global y-derivatives at 
 * the interfaces.
 */
template <
        class PatchSeqMin,
        class PatchSeqMax,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type,
        class... Patches>
void check_all_y_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_range_slices_dy,
        std::tuple<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>...>
                coord_transforms)
{
    (check_y_derivatives<Patches, PatchSeqMin, PatchSeqMax>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             typename Patches::IdxRange2(idx_ranges.template get<Patches>()),
             idx_range_slices_dy.template get<Patches>(),
             std::get<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>>(
                     coord_transforms)),
     ...);
}


/** 
 * @brief Check agreement between the computed cross-derivatives and the global cross-derivatives at 
 * the interfaces for a given patch.
 */
template <
        class Patch,
        class PatchSeqMin,
        class PatchSeqMax,
        class CoordTransformType,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2>
void check_xy_derivatives(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        typename Patch::IdxRange12 const& idx_range,
        IdxRange1SliceOnPatch<Patch> const& idx_range_slice_dx,
        IdxRange2SliceOnPatch<Patch> const& idx_range_slice_dy,
        CoordTransformType const& coord_transform = CoordTransformType())
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    const int sign_x = !coord_transform.m_is_reverse_x - coord_transform.m_is_reverse_x;
    const int sign_y = !coord_transform.m_is_reverse_y - coord_transform.m_is_reverse_y;
    const int sign_xy = sign_x * sign_y;

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
            coord_transform.get_global_coord(ddc::coordinate(idx_min_min)));
    Coord<Xg, Yg> interface_coord_max_min(
            coord_transform.get_global_coord(ddc::coordinate(idx_max_min)));
    Coord<Xg, Yg> interface_coord_min_max(
            coord_transform.get_global_coord(ddc::coordinate(idx_min_max)));
    Coord<Xg, Yg> interface_coord_max_max(
            coord_transform.get_global_coord(ddc::coordinate(idx_max_max)));

    double global_deriv_min_min
            = sign_xy * evaluator_g.deriv_1_and_2(interface_coord_min_min, function_g_coef);
    double global_deriv_max_min
            = sign_xy * evaluator_g.deriv_1_and_2(interface_coord_max_min, function_g_coef);
    double global_deriv_min_max
            = sign_xy * evaluator_g.deriv_1_and_2(interface_coord_min_max, function_g_coef);
    double global_deriv_max_max
            = sign_xy * evaluator_g.deriv_1_and_2(interface_coord_max_max, function_g_coef);

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

/** 
 * @brief Check agreement between the computed cross-derivatives and the global cross-derivatives at 
 * the interfaces.
 */
template <
        class PatchSeqMin,
        class PatchSeqMax,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type,
        class... Patches>
void check_all_xy_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_range_slices_dx,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_range_slices_dy,
        std::tuple<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>...>
                coord_transforms)
{
    (check_xy_derivatives<Patches, PatchSeqMin, PatchSeqMax>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             idx_ranges.template get<Patches>(),
             idx_range_slices_dx.template get<Patches>(),
             idx_range_slices_dy.template get<Patches>(),
             std::get<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>>(
                     coord_transforms)),
     ...);
}


/** 
 * @brief Check agreement between the local splines defined with the computed derivatives 
 * and the global spline for a given patch.  
 */
template <
        class Patch,
        class PatchSeqMin,
        class PatchSeqMax,
        class CoordTransformType,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2>
void check_spline_representation_agreement(
        typename Patch::IdxRange12 const& idx_range_xy,
        IdxRange1SliceOnPatch<Patch> const& idx_range_slice_x,
        IdxRange2SliceOnPatch<Patch> const& idx_range_slice_y,
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        CoordTransformType const& coord_transform = CoordTransformType())
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

    DerivFieldSplineBuilder<
            HostExecSpace,
            typename Patch::BSplines1,
            typename Patch::BSplines2,
            typename Patch::Grid1,
            typename Patch::Grid2,
            BoundCondXmin,
            BoundCondXmax,
            BoundCondYmin,
            BoundCondYmax>
            apply_builder(builder);

    SplineCoeffMemOnPatch_2D_host<Patch> function_coef_alloc(
            builder.batched_spline_domain(idx_range_xy));
    SplineCoeffOnPatch_2D_host<Patch> function_coef = get_field(function_coef_alloc);

    apply_builder(function_coef, function_and_derivs);

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
        Coord<Xg, Yg> eval_point_g(coord_transform.get_global_coord(eval_points(idx)));
        double local_spline = evaluator(eval_points(idx), get_const_field(function_coef));
        double global_spline = evaluator_g(eval_point_g, get_const_field(function_g_coef));
        EXPECT_NEAR(local_spline, global_spline, 1e-14);
    });
}

/** 
 * @brief Check agreement between the local splines defined with the computed derivatives 
 * and the global spline. 
 */
template <
        class PatchSeqMin,
        class PatchSeqMax,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type,
        class... Patches>
void check_all_spline_representation_agreement(
        MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
        MultipatchType<IdxRange1SliceOnPatch, Patches...> const& idx_range_slices_x,
        MultipatchType<IdxRange2SliceOnPatch, Patches...> const& idx_range_slices_y,
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        std::tuple<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>...>
                coord_transforms)
{
    (check_spline_representation_agreement<Patches, PatchSeqMin, PatchSeqMax>(
             idx_ranges.template get<Patches>(),
             idx_range_slices_x.template get<Patches>(),
             idx_range_slices_y.template get<Patches>(),
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             get_const_field(function_g_coef),
             std::get<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>>(
                     coord_transforms)),
     ...);
};
