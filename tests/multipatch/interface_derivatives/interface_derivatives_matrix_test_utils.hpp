
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "multipatch_field.hpp"
#include "spline_builder_deriv_field_2d.hpp"
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

    const bool m_is_x_loc_well_oriented;
    const bool m_is_y_loc_well_oriented;
    const bool are_x_y_loc_x_y_glob;

    const Coord<X_loc> m_x_min;
    const Coord<X_loc> m_x_max;
    const Coord<Y_loc> m_y_min;
    const Coord<Y_loc> m_y_max;

    /**
     * @brief Instantiate the coordinate transformator.
     * By default, the transformation is identity.
     */
    explicit CoordTransform(
            bool is_x_loc_well_oriented = true,
            bool is_y_loc_well_oriented = true,
            bool are_exchange_x_y = true,
            Coord<X_loc> x_min = Coord<X_loc>(0),
            Coord<X_loc> x_max = Coord<X_loc>(1),
            Coord<Y_loc> y_min = Coord<Y_loc>(0),
            Coord<Y_loc> y_max = Coord<Y_loc>(1))
        : m_is_x_loc_well_oriented(is_x_loc_well_oriented)
        , m_is_y_loc_well_oriented(is_y_loc_well_oriented)
        , are_x_y_loc_x_y_glob(are_exchange_x_y)
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

        if (!m_is_x_loc_well_oriented) {
            x_loc = m_x_min + m_x_max - x_loc;
        }
        if (!m_is_y_loc_well_oriented) {
            y_loc = m_y_min + m_y_max - y_loc;
        }
        if (!are_x_y_loc_x_y_glob) {
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
    ddc::host_for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
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
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef,
        CoordTransformType const& coord_transform = CoordTransformType())
{
    using DerivX = ddc::Deriv<typename Patch::Dim1>;
    using DerivY = ddc::Deriv<typename Patch::Dim2>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    IdxRange<GridX> idx_range_x(get_idx_range(function_and_derivs));
    IdxRange<GridY> idx_range_y(get_idx_range(function_and_derivs));

    Idx<DerivY> idx_dy(Idx<DerivY>(1));
    Idx<GridY> idx_ymin(idx_range_y.front());
    Idx<GridY> idx_ymax(idx_range_y.back());

    int const sign
            = coord_transform.m_is_y_loc_well_oriented - !coord_transform.m_is_y_loc_well_oriented;

    ddc::host_for_each(idx_range_x, [&](Idx<GridX> const& idx_par) {
        Idx<GridX, GridY> idx_min(idx_par, idx_ymin);
        Idx<GridX, GridY> idx_max(idx_par, idx_ymax);
        Coord<Xg, Yg> interface_coord_min(
                coord_transform.get_global_coord(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(
                coord_transform.get_global_coord(ddc::coordinate(idx_max)));

        if (coord_transform.are_x_y_loc_x_y_glob) {
            Idx<ddc::Deriv<Yg>> idx_deriv(1);
            function_and_derivs(idx_dy, idx_ymin, idx_par)
                    = sign
                      * evaluator_g.deriv(idx_deriv, interface_coord_min, const_function_g_coef);
            function_and_derivs(idx_dy, idx_ymax, idx_par)
                    = sign
                      * evaluator_g.deriv(idx_deriv, interface_coord_max, const_function_g_coef);
        } else {
            Idx<ddc::Deriv<Xg>> idx_deriv(1);
            function_and_derivs(idx_dy, idx_ymin, idx_par)
                    = sign
                      * evaluator_g.deriv(idx_deriv, interface_coord_min, const_function_g_coef);
            function_and_derivs(idx_dy, idx_ymax, idx_par)
                    = sign
                      * evaluator_g.deriv(idx_deriv, interface_coord_max, const_function_g_coef);
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
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef,
        CoordTransformType const& coord_transform = CoordTransformType())
{
    using DerivX = ddc::Deriv<typename Patch::Dim1>;
    using DerivY = ddc::Deriv<typename Patch::Dim2>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    IdxRange<GridX> idx_range_x(get_idx_range(function_and_derivs));
    IdxRange<GridY> idx_range_y(get_idx_range(function_and_derivs));

    Idx<DerivX> idx_dx(Idx<DerivX>(1));
    Idx<GridX> idx_xmin(idx_range_x.front());
    Idx<GridX> idx_xmax(idx_range_x.back());

    int const sign
            = coord_transform.m_is_x_loc_well_oriented - !coord_transform.m_is_x_loc_well_oriented;

    ddc::host_for_each(idx_range_y, [&](Idx<GridY> const& idx_par) {
        Idx<GridX, GridY> idx_min(idx_xmin, idx_par);
        Idx<GridX, GridY> idx_max(idx_xmax, idx_par);
        Coord<Xg, Yg> interface_coord_min(
                coord_transform.get_global_coord(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(
                coord_transform.get_global_coord(ddc::coordinate(idx_max)));

        if (coord_transform.are_x_y_loc_x_y_glob) {
            Idx<ddc::Deriv<Xg>> idx_deriv(1);
            function_and_derivs(idx_dx, idx_xmin, idx_par)
                    = sign
                      * evaluator_g.deriv(idx_deriv, interface_coord_min, const_function_g_coef);
            function_and_derivs(idx_dx, idx_xmax, idx_par)
                    = sign
                      * evaluator_g.deriv(idx_deriv, interface_coord_max, const_function_g_coef);
        } else {
            Idx<ddc::Deriv<Yg>> idx_deriv(1);
            function_and_derivs(idx_dx, idx_xmin, idx_par)
                    = sign
                      * evaluator_g.deriv(idx_deriv, interface_coord_min, const_function_g_coef);
            function_and_derivs(idx_dx, idx_xmax, idx_par)
                    = sign
                      * evaluator_g.deriv(idx_deriv, interface_coord_max, const_function_g_coef);
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

    IdxRange<GridX> idx_range_x(get_idx_range(function_and_derivs));
    IdxRange<GridY> idx_range_y(get_idx_range(function_and_derivs));

    Idx<DerivX> idx_dx(Idx<DerivX>(1));
    Idx<DerivY> idx_dy(Idx<DerivY>(1));
    Idx<GridX> idx_xmin(idx_range_x.front());
    Idx<GridX> idx_xmax(idx_range_x.back());
    Idx<GridY> idx_ymin(idx_range_y.front());
    Idx<GridY> idx_ymax(idx_range_y.back());

    IdxdXXdYY idx_cross_deriv_min_min(idx_dx, idx_xmin, idx_dy, idx_ymin);
    IdxdXXdYY idx_cross_deriv_max_min(idx_dx, idx_xmax, idx_dy, idx_ymin);
    IdxdXXdYY idx_cross_deriv_min_max(idx_dx, idx_xmin, idx_dy, idx_ymax);
    IdxdXXdYY idx_cross_deriv_max_max(idx_dx, idx_xmax, idx_dy, idx_ymax);

    Coord<Xg, Yg> coord_min_min = coord_transform.get_global_coord(
            typename Patch::Coord12(coord_transform.m_x_min, coord_transform.m_y_min));
    Coord<Xg, Yg> coord_max_min = coord_transform.get_global_coord(
            typename Patch::Coord12(coord_transform.m_x_max, coord_transform.m_y_min));
    Coord<Xg, Yg> coord_min_max = coord_transform.get_global_coord(
            typename Patch::Coord12(coord_transform.m_x_min, coord_transform.m_y_max));
    Coord<Xg, Yg> coord_max_max = coord_transform.get_global_coord(
            typename Patch::Coord12(coord_transform.m_x_max, coord_transform.m_y_max));

    int const sign_x
            = coord_transform.m_is_x_loc_well_oriented - !coord_transform.m_is_x_loc_well_oriented;
    int const sign_y
            = coord_transform.m_is_y_loc_well_oriented - !coord_transform.m_is_y_loc_well_oriented;
    int const sign_xy = sign_x * sign_y;

    Idx<ddc::Deriv<Xg>, ddc::Deriv<Yg>> idx_deriv(1, 1);
    function_and_derivs(idx_cross_deriv_min_min)
            = sign_xy * evaluator_g.deriv(idx_deriv, coord_min_min, const_function_g_coef);
    function_and_derivs(idx_cross_deriv_max_min)
            = sign_xy * evaluator_g.deriv(idx_deriv, coord_max_min, const_function_g_coef);
    function_and_derivs(idx_cross_deriv_min_max)
            = sign_xy * evaluator_g.deriv(idx_deriv, coord_min_max, const_function_g_coef);
    function_and_derivs(idx_cross_deriv_max_max)
            = sign_xy * evaluator_g.deriv(idx_deriv, coord_max_max, const_function_g_coef);
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
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef,
        std::tuple<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>...>
                coord_transforms)
{
    using PatchSeq = ddc::detail::TypeSeq<Patches...>;
    (initialise_cross_derivatives<Patches>(
             functions_and_derivs.template get<Patches>(),
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
    ddc::host_for_each(idx_range, [&](typename Patch::Idx12 const& idx) {
        IdxStep<Grid1> idx_x = coord_transform.m_is_x_loc_well_oriented
                                       ? Idx<Grid1>(idx) - IdxRange<Grid1>(idx_range).front()
                                       : IdxRange<Grid1>(idx_range).back() - Idx<Grid1>(idx);
        IdxStep<Grid2> idx_y = coord_transform.m_is_y_loc_well_oriented
                                       ? Idx<Grid2>(idx) - IdxRange<Grid2>(idx_range).front()
                                       : IdxRange<Grid2>(idx_range).back() - Idx<Grid2>(idx);
        Idx<GridXg, GridYg> idx_g
                = coord_transform.are_x_y_loc_x_y_glob
                          ? Idx<GridXg, GridYg>(idx_x.value() + x_shift, idx_y.value() + y_shift)
                          : Idx<GridXg, GridYg>(idx_y.value() + x_shift, idx_x.value() + y_shift);

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
        CoordTransformType const& coord_transform = CoordTransformType(),
        double const TOL = 5e-14)
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    const int sign
            = coord_transform.m_is_x_loc_well_oriented - !coord_transform.m_is_x_loc_well_oriented;

    IdxRange<GridX> idx_range_x(get_idx_range(function_and_derivs));
    IdxRange<GridY> idx_range_y(get_idx_range(function_and_derivs));

    Idx<DerivX> idx_deriv(1);
    Idx<GridX> idx_xmin(idx_range_x.front());
    Idx<GridX> idx_xmax(idx_range_x.back());

    ddc::host_for_each(idx_range_y, [&](typename Patch::Idx2 const& idx_y) {
        typename Patch::Idx12 idx_min(idx_xmin, idx_y);
        typename Patch::Idx12 idx_max(idx_xmax, idx_y);
        Coord<Xg, Yg> interface_coord_min(
                coord_transform.get_global_coord(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(
                coord_transform.get_global_coord(ddc::coordinate(idx_max)));

        double global_deriv_min;
        double global_deriv_max;
        if (coord_transform.are_x_y_loc_x_y_glob) {
            Idx<ddc::Deriv<Xg>> idx_deriv(1);
            global_deriv_min
                    = sign * evaluator_g.deriv(idx_deriv, interface_coord_min, function_g_coef);
            global_deriv_max
                    = sign * evaluator_g.deriv(idx_deriv, interface_coord_max, function_g_coef);
        } else {
            Idx<ddc::Deriv<Yg>> idx_deriv(1);
            global_deriv_min
                    = sign * evaluator_g.deriv(idx_deriv, interface_coord_min, function_g_coef);
            global_deriv_max
                    = sign * evaluator_g.deriv(idx_deriv, interface_coord_max, function_g_coef);
        }
        EXPECT_NEAR(function_and_derivs(idx_deriv, idx_xmin, idx_y), global_deriv_min, TOL);
        EXPECT_NEAR(function_and_derivs(idx_deriv, idx_xmax, idx_y), global_deriv_max, TOL);
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
        std::tuple<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>...>
                coord_transforms,
        double TOL = 5e-14)
{
    (check_x_derivatives<Patches>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             std::get<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>>(
                     coord_transforms),
             TOL),
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
        CoordTransformType const& coord_transform = CoordTransformType(),
        double const TOL = 5e-14)
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    using DerivY = typename ddc::Deriv<typename Patch::Dim2>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    const int sign
            = coord_transform.m_is_y_loc_well_oriented - !coord_transform.m_is_y_loc_well_oriented;

    IdxRange<GridX> idx_range_x(get_idx_range(function_and_derivs));
    IdxRange<GridY> idx_range_y(get_idx_range(function_and_derivs));

    Idx<DerivY> idx_deriv(1);
    Idx<GridY> idx_ymin = idx_range_y.front();
    Idx<GridY> idx_ymax = idx_range_y.back();

    ddc::host_for_each(idx_range_x, [&](typename Patch::Idx1 const& idx_x) {
        typename Patch::Idx12 idx_min(idx_x, idx_ymin);
        typename Patch::Idx12 idx_max(idx_x, idx_ymax);
        Coord<Xg, Yg> interface_coord_min(
                coord_transform.get_global_coord(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(
                coord_transform.get_global_coord(ddc::coordinate(idx_max)));

        double global_deriv_min;
        double global_deriv_max;
        if (coord_transform.are_x_y_loc_x_y_glob) {
            Idx<ddc::Deriv<Yg>> idx_deriv(1);
            global_deriv_min
                    = sign * evaluator_g.deriv(idx_deriv, interface_coord_min, function_g_coef);
            global_deriv_max
                    = sign * evaluator_g.deriv(idx_deriv, interface_coord_max, function_g_coef);
        } else {
            Idx<ddc::Deriv<Xg>> idx_deriv(1);
            global_deriv_min
                    = sign * evaluator_g.deriv(idx_deriv, interface_coord_min, function_g_coef);
            global_deriv_max
                    = sign * evaluator_g.deriv(idx_deriv, interface_coord_max, function_g_coef);
        }

        // For Patches in PatchSeqMin, we defined ddc::BoundCond::GREVILLE the local lower Y-boundary,
        // we don't need the y-derivatives for y = ymin. Their values are not checked.
        if constexpr (!ddc::in_tags_v<Patch, PatchSeqMin>) {
            EXPECT_NEAR(function_and_derivs(idx_deriv, idx_ymin, idx_x), global_deriv_min, TOL);
        }
        // For Patches in PatchSeqMax, we defined ddc::BoundCond::GREVILLE the local lower Y-boundary,
        // we don't need the y-derivatives for y = ymax. Their values are not checked.
        if constexpr (!ddc::in_tags_v<Patch, PatchSeqMax>) {
            EXPECT_NEAR(function_and_derivs(idx_deriv, idx_ymax, idx_x), global_deriv_max, TOL);
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
        std::tuple<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>...>
                coord_transforms,
        double const TOL = 5e-14)
{
    (check_y_derivatives<Patches, PatchSeqMin, PatchSeqMax>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             std::get<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>>(
                     coord_transforms),
             TOL),
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
        CoordTransformType const& coord_transform = CoordTransformType(),
        double const TOL = 5e-13)
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    const int sign_x
            = coord_transform.m_is_x_loc_well_oriented - !coord_transform.m_is_x_loc_well_oriented;
    const int sign_y
            = coord_transform.m_is_y_loc_well_oriented - !coord_transform.m_is_y_loc_well_oriented;
    const int sign_xy = sign_x * sign_y;

    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;
    using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
    using DerivY = typename ddc::Deriv<typename Patch::Dim2>;
    using IdxdXXdYY = Idx<DerivX, GridX, DerivY, GridY>;

    IdxRange<GridX> idx_range_x(get_idx_range(function_and_derivs));
    IdxRange<GridY> idx_range_y(get_idx_range(function_and_derivs));

    Idx<DerivX> idx_dx(1);
    Idx<DerivY> idx_dy(1);

    Idx<GridX> idx_xmin(idx_range_x.front());
    Idx<GridX> idx_xmax(idx_range_x.back());
    Idx<GridY> idx_ymin(idx_range_y.front());
    Idx<GridY> idx_ymax(idx_range_y.back());

    IdxdXXdYY idx_cross_deriv_min_min(idx_dx, idx_xmin, idx_dy, idx_ymin);
    IdxdXXdYY idx_cross_deriv_max_min(idx_dx, idx_xmax, idx_dy, idx_ymin);
    IdxdXXdYY idx_cross_deriv_min_max(idx_dx, idx_xmin, idx_dy, idx_ymax);
    IdxdXXdYY idx_cross_deriv_max_max(idx_dx, idx_xmax, idx_dy, idx_ymax);

    Idx<GridX, GridY> idx_min_min(idx_xmin, idx_ymin);
    Idx<GridX, GridY> idx_max_min(idx_xmax, idx_ymin);
    Idx<GridX, GridY> idx_min_max(idx_xmin, idx_ymax);
    Idx<GridX, GridY> idx_max_max(idx_xmax, idx_ymax);

    Coord<Xg, Yg> interface_coord_min_min(
            coord_transform.get_global_coord(ddc::coordinate(idx_min_min)));
    Coord<Xg, Yg> interface_coord_max_min(
            coord_transform.get_global_coord(ddc::coordinate(idx_max_min)));
    Coord<Xg, Yg> interface_coord_min_max(
            coord_transform.get_global_coord(ddc::coordinate(idx_min_max)));
    Coord<Xg, Yg> interface_coord_max_max(
            coord_transform.get_global_coord(ddc::coordinate(idx_max_max)));

    Idx<ddc::Deriv<Xg>, ddc::Deriv<Yg>> idx_deriv(1, 1);
    double global_deriv_min_min
            = sign_xy * evaluator_g.deriv(idx_deriv, interface_coord_min_min, function_g_coef);
    double global_deriv_max_min
            = sign_xy * evaluator_g.deriv(idx_deriv, interface_coord_max_min, function_g_coef);
    double global_deriv_min_max
            = sign_xy * evaluator_g.deriv(idx_deriv, interface_coord_min_max, function_g_coef);
    double global_deriv_max_max
            = sign_xy * evaluator_g.deriv(idx_deriv, interface_coord_max_max, function_g_coef);

    // For Patches in PatchSeqMin, we defined ddc::BoundCond::GREVILLE the local lower Y-boundary,
    // we don't need the cross-derivatives for y = ymin. Their value is not checked.
    if constexpr (!ddc::in_tags_v<Patch, PatchSeqMin>) {
        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_min_min), global_deriv_min_min, TOL);
        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_max_min), global_deriv_max_min, TOL);
    }
    // For Patches in PatchSeqMax, we defined ddc::BoundCond::GREVILLE the local upper Y-boundary,
    // we don't need the cross-derivatives for y = ymax. Their value is not checked.
    if constexpr (!ddc::in_tags_v<Patch, PatchSeqMax>) {
        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_min_max), global_deriv_min_max, TOL);
        EXPECT_NEAR(function_and_derivs(idx_cross_deriv_max_max), global_deriv_max_max, TOL);
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
        std::tuple<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>...>
                coord_transforms,
        double const TOL = 5e-13)
{
    (check_xy_derivatives<Patches, PatchSeqMin, PatchSeqMax>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             std::get<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>>(
                     coord_transforms),
             TOL),
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
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        CoordTransformType const& coord_transform = CoordTransformType(),
        double const TOL = 1e-14)
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

    IdxRange<GridX, GridY> idx_range_xy = get_idx_range(function_and_derivs);

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

    SplineBuliderDerivField2D<
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
    ddc::host_for_each(idx_range_xy, [&](typename Patch::Idx12 const idx) {
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
    ddc::host_for_each(idx_range_xy, [&](typename Patch::Idx12 const idx) {
        Coord<Xg, Yg> eval_point_g(coord_transform.get_global_coord(eval_points(idx)));
        double local_spline = evaluator(eval_points(idx), get_const_field(function_coef));
        double global_spline = evaluator_g(eval_point_g, get_const_field(function_g_coef));
        EXPECT_NEAR(local_spline, global_spline, TOL);
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
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        std::tuple<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>...>
                coord_transforms,
        double const TOL = 5e-14)
{
    (check_spline_representation_agreement<Patches, PatchSeqMin, PatchSeqMax>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             get_const_field(function_g_coef),
             std::get<CoordTransform<Xg, Yg, typename Patches::Dim1, typename Patches::Dim2>>(
                     coord_transforms),
             TOL),
     ...);
};
