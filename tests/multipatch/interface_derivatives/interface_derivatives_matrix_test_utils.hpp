
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
template <bool IsXTransform, class CoordTransformType>
auto get_1d_transform(CoordTransformType transform)
{
    if constexpr (IsXTransform) {
        return transform.x_transform;
    } else {
        return transform.y_transform;
    }
}

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
template <class Xg, class Yg, class Grid1, class Grid2, class Layout>
void initialise_2D_function(DField<IdxRange<Grid1, Grid2>, Kokkos::HostSpace, Layout> function)
{
    ddc::host_for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        Coord<Xg, Yg> equiv_global_coord = ddc::coordinate(idx);
        double const xg = ddc::get<Xg>(equiv_global_coord);
        double const yg = ddc::get<Yg>(equiv_global_coord);
        function(idx) = Kokkos::cos(xg * 2. / 3. * M_PI + 0.25) * Kokkos::sin(yg);
    });
}

/// @brief Initialise the function with f(x,y) = cos(2/3*pi*x)sin(y).
template <class Xg, class Yg, class Grid1, class Grid2, class Layout, class CoordTransformType>
void initialise_2D_function(
        DField<IdxRange<Grid1, Grid2>, Kokkos::HostSpace, Layout> function,
        CoordTransformType const& coord_transform)
{
    auto coord_transform_l_to_g = coord_transform.coord_transform.get_inverse_mapping();
    ddc::host_for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        Coord<Xg, Yg> equiv_global_coord(coord_transform_l_to_g(ddc::coordinate(idx)));
        double const xg = ddc::get<Xg>(equiv_global_coord);
        double const yg = ddc::get<Yg>(equiv_global_coord);
        function(idx) = Kokkos::cos(xg * 2. / 3. * M_PI + 0.25) * Kokkos::sin(yg);
    });
}

template <class Xg, class Yg, class... Patches, class... CoordTransformType, size_t... I>
void initialise_all_functions(
        MultipatchField<DerivFieldOnPatch_host, Patches...> const& functions_and_derivs,
        std::tuple<CoordTransformType...> coord_transforms)
{
    static_assert(sizeof...(Patches) == sizeof...(CoordTransformType));
    using PatchSeq = ddc::detail::TypeSeq<Patches...>;
    (initialise_2D_function<Xg, Yg, typename Patches::Grid1, typename Patches::Grid2>(
             (functions_and_derivs.template get<Patches>()).get_values_field(),
             std::get<ddc::type_seq_rank_v<Patches, PatchSeq>>(coord_transforms)),
     ...);
}

template <
        class Patch,
        class GridDeriv,
        class CoordTransformType,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class Xg = typename BSplinesXg::continuous_dimension_type,
        class Yg = typename BSplinesYg::continuous_dimension_type>
void initialise_derivatives(
        Idx<GridDeriv> idx_deriv_pos,
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef,
        CoordTransformType coord_transform)
{
    using DerivDimG = typename GridDeriv::continuous_dimension_type;
    static_assert(ddc::in_tags_v<DerivDimG, ddc::detail::TypeSeq<Xg, Yg>>);

    auto coord_transform_g_to_l = get_1d_transform<std::is_same_v<DerivDimG, Xg>>(coord_transform);
    auto coord_transform_l_to_g = coord_transform.coord_transform.get_inverse_mapping();

    using CoordAlongDeriv = decltype(coord_transform_g_to_l)::CoordResult;
    using DimAlongDeriv = ddc::type_seq_element_t<0, ddc::to_type_seq_t<CoordAlongDeriv>>;
    using GridAlongDeriv = std::conditional_t<
            std::is_same_v<typename Patch::Grid1::continuous_dimension_type, DimAlongDeriv>,
            typename Patch::Grid1,
            typename Patch::Grid2>;
    using GridPerpToDeriv = std::conditional_t<
            std::is_same_v<typename Patch::Grid1, GridAlongDeriv>,
            typename Patch::Grid2,
            typename Patch::Grid1>;

    Coord<DerivDimG> deriv_coord_glob = ddc::coordinate(idx_deriv_pos);
    Coord<DimAlongDeriv> deriv_coord = coord_transform_g_to_l(deriv_coord_glob);

    IdxRange<GridAlongDeriv> idx_range_derivs(get_idx_range(function_and_derivs));
    IdxRange<GridPerpToDeriv> idx_range_perp(get_idx_range(function_and_derivs));

    Idx<ddc::Deriv<DerivDimG>> idx_deriv(2);
    Idx<GridAlongDeriv> idx_deriv_pos_local
            = (std::fabs(deriv_coord - ddc::coordinate(idx_range_derivs.front())) < 1e-14)
                      ? idx_range_derivs.front()
                      : idx_range_derivs.back();

    ddc::host_for_each(idx_range_perp, [&](Idx<GridPerpToDeriv> const& idx_perp) {
        Coord<Xg, Yg> coord_glob(deriv_coord_glob, coord_transform_l_to_g(ddc::coordinate(idx_perp)));
        function_and_derivs(idx_deriv, idx_deriv_pos_local, idx_perp)
                = coord_transform_g_to_l.jacobian(deriv_coord_glob)
                  * evaluator_g.deriv(idx_deriv, coord_glob, const_function_g_coef);
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
        CoordTransformType const& coord_transform_handler)
{
    using DerivX = ddc::Deriv<typename Patch::Dim1>;
    using DerivY = ddc::Deriv<typename Patch::Dim2>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    IdxRange<GridX> idx_range_x(get_idx_range(function_and_derivs));
    IdxRange<GridY> idx_range_y(get_idx_range(function_and_derivs));

    Idx<DerivX> idx_dx(Idx<DerivX>(1));
    Idx<DerivY> idx_dy(Idx<DerivY>(1));
    Idx<GridX> idx_xmin(idx_range_x.front());
    Idx<GridX> idx_xmax(idx_range_x.back());
    Idx<GridY> idx_ymin(idx_range_y.front());
    Idx<GridY> idx_ymax(idx_range_y.back());

    typename Patch::Coord12 coord_min_min(ddc::coordinate(idx_xmin), ddc::coordinate(idx_ymin));
    typename Patch::Coord12 coord_max_min(ddc::coordinate(idx_xmax), ddc::coordinate(idx_ymin));
    typename Patch::Coord12 coord_min_max(ddc::coordinate(idx_xmin), ddc::coordinate(idx_ymax));
    typename Patch::Coord12 coord_max_max(ddc::coordinate(idx_xmax), ddc::coordinate(idx_ymax));

    auto coord_transform = coord_transform_handler.coord_transform;
    auto coord_transform_l_to_g = coord_transform_handler.coord_transform.get_inverse_mapping();

    Coord<Xg, Yg> coord_min_min_g(coord_transform_l_to_g(coord_min_min));
    Coord<Xg, Yg> coord_max_min_g(coord_transform_l_to_g(coord_max_min));
    Coord<Xg, Yg> coord_min_max_g(coord_transform_l_to_g(coord_min_max));
    Coord<Xg, Yg> coord_max_max_g(coord_transform_l_to_g(coord_max_max));

    Idx<ddc::Deriv<Xg>, ddc::Deriv<Yg>> idx_deriv(1, 1);
    function_and_derivs(idx_dx, idx_xmin, idx_dy, idx_ymin)
            = coord_transform.jacobian(coord_min_min_g)
              * evaluator_g.deriv(idx_deriv, coord_min_min_g, const_function_g_coef);
    function_and_derivs(idx_dx, idx_xmax, idx_dy, idx_ymin)
            = coord_transform.jacobian(coord_max_min_g)
              * evaluator_g.deriv(idx_deriv, coord_max_min_g, const_function_g_coef);
    function_and_derivs(idx_dx, idx_xmin, idx_dy, idx_ymax)
            = coord_transform.jacobian(coord_min_max_g)
              * evaluator_g.deriv(idx_deriv, coord_min_max_g, const_function_g_coef);
    function_and_derivs(idx_dx, idx_xmax, idx_dy, idx_ymax)
            = coord_transform.jacobian(coord_max_max_g)
              * evaluator_g.deriv(idx_deriv, coord_max_max_g, const_function_g_coef);
}

/// @brief Initialise all the cross-derivatives of the given DerivFields from the global spline.
template <
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2,
        class... Patches,
        class... CoordTransformType>
void initialise_all_cross_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& const_function_g_coef,
        std::tuple<CoordTransformType...> coord_transforms)
{
    using PatchSeq = ddc::detail::TypeSeq<Patches...>;
    (initialise_cross_derivatives<Patches>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             const_function_g_coef,
             ddc::type_seq_rank_v<Patches, PatchSeq>,
             std::get<ddc::type_seq_rank_v<Patches, PatchSeq>>(coord_transforms)),
     ...);
}

// -----------------------------------------------------------------------------------------------
// CHECK OPERATORS -------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------


/** 
 * @brief Check agreement between the computed x-derivatives and the global x-derivatives at 
 * the interfaces for a given patch.
 */
template <
        class DerivDirection,
        class Patch,
        class CoordTransformType,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2>
void check_derivatives(
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        CoordTransformType const& coord_transform,
        double const TOL = 5e-14)
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    using DerivG = ddc::Deriv<std::conditional_t<
            std::is_same_v<typename CoordTransformType::XTransform::CoordResult, DerivX>,
            Xg,
            Yg>>;

    IdxRange<GridX> idx_range_x(get_idx_range(function_and_derivs));
    IdxRange<GridY> idx_range_y(get_idx_range(function_and_derivs));

    Idx<DerivX> idx_deriv(1);
    Idx<GridX> idx_xmin(idx_range_x.front());
    Idx<GridX> idx_xmax(idx_range_x.back());

    auto coord_transform_g_to_l = coord_transform.coord_transform;
    auto coord_transform_l_to_g = coord_transform.coord_transform.get_inverse_mapping();

    ddc::host_for_each(idx_range_y, [&](typename Patch::Idx2 const& idx_y) {
        typename Patch::Idx12 idx_min(idx_xmin, idx_y);
        typename Patch::Idx12 idx_max(idx_xmax, idx_y);
        Coord<Xg, Yg> interface_coord_min(coord_transform_l_to_g(ddc::coordinate(idx_min)));
        Coord<Xg, Yg> interface_coord_max(coord_transform_l_to_g(ddc::coordinate(idx_max)));

        double global_deriv_min;
        double global_deriv_max;
        Idx<DerivG> idx_deriv(1);
        global_deriv_min = coord_transform_g_to_l.jacobian(interface_coord_min)
                           * evaluator_g.deriv(idx_deriv, interface_coord_min, function_g_coef);
        global_deriv_max = coord_transform_g_to_l.jacobian(interface_coord_max)
                           * evaluator_g.deriv(idx_deriv, interface_coord_max, function_g_coef);
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
        class... Patches,
        class... CoordTransformType>
void check_all_x_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        std::tuple<CoordTransformType...> coord_transforms,
        double TOL = 5e-14)
{
    using PatchSeq = ddc::detail::TypeSeq<Patches...>;
    (check_derivatives<typename Patches::Dim1, Patches>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             std::get<ddc::type_seq_rank_v<Patches, PatchSeq>>(coord_transforms),
             TOL),
     ...);
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
        class... Patches,
        class... CoordTransformType>
void check_all_y_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        std::tuple<CoordTransformType...> coord_transforms,
        double const TOL = 5e-14)
{
    using PatchSeq = ddc::detail::TypeSeq<Patches...>;
    (check_derivatives<typename Patches::Dim2, Patches>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             std::get<ddc::type_seq_rank_v<Patches, PatchSeq>>(coord_transforms),
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
        CoordTransformType const& coord_transform,
        double const TOL = 5e-13)
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;
    using DerivX = typename ddc::Deriv<typename Patch::Dim1>;
    using DerivY = typename ddc::Deriv<typename Patch::Dim2>;

    IdxRange<GridX> idx_range_x(get_idx_range(function_and_derivs));
    IdxRange<GridY> idx_range_y(get_idx_range(function_and_derivs));

    Idx<DerivX> idx_dx(1);
    Idx<DerivY> idx_dy(1);

    Idx<GridX> idx_xmin(idx_range_x.front());
    Idx<GridX> idx_xmax(idx_range_x.back());
    Idx<GridY> idx_ymin(idx_range_y.front());
    Idx<GridY> idx_ymax(idx_range_y.back());

    Idx<GridX, GridY> idx_min_min(idx_xmin, idx_ymin);
    Idx<GridX, GridY> idx_max_min(idx_xmax, idx_ymin);
    Idx<GridX, GridY> idx_min_max(idx_xmin, idx_ymax);
    Idx<GridX, GridY> idx_max_max(idx_xmax, idx_ymax);

    auto coord_transform_g_to_l = coord_transform.coord_transform;
    auto coord_transform_l_to_g = coord_transform.coord_transform.get_inverse_mapping();

    typename Patch::Coord12 coord_min_min(ddc::coordinate(idx_min_min));
    typename Patch::Coord12 coord_max_min(ddc::coordinate(idx_max_min));
    typename Patch::Coord12 coord_min_max(ddc::coordinate(idx_min_max));
    typename Patch::Coord12 coord_max_max(ddc::coordinate(idx_max_max));

    Coord<Xg, Yg> coord_min_min_g(coord_transform_l_to_g(coord_min_min));
    Coord<Xg, Yg> coord_max_min_g(coord_transform_l_to_g(coord_max_min));
    Coord<Xg, Yg> coord_min_max_g(coord_transform_l_to_g(coord_min_max));
    Coord<Xg, Yg> coord_max_max_g(coord_transform_l_to_g(coord_max_max));

    Idx<ddc::Deriv<Xg>, ddc::Deriv<Yg>> idx_deriv(1, 1);
    double global_deriv_min_min = coord_transform_g_to_l.jacobian(coord_min_min_g)
                                  * evaluator_g.deriv(idx_deriv, coord_min_min_g, function_g_coef);
    double global_deriv_max_min = coord_transform_g_to_l.jacobian(coord_max_min_g)
                                  * evaluator_g.deriv(idx_deriv, coord_max_min_g, function_g_coef);
    double global_deriv_min_max = coord_transform_g_to_l.jacobian(coord_min_max_g)
                                  * evaluator_g.deriv(idx_deriv, coord_min_max_g, function_g_coef);
    double global_deriv_max_max = coord_transform_g_to_l.jacobian(coord_max_max_g)
                                  * evaluator_g.deriv(idx_deriv, coord_max_max_g, function_g_coef);

    // For Patches in PatchSeqMin, we defined ddc::BoundCond::GREVILLE the local lower Y-boundary,
    // we don't need the cross-derivatives for y = ymin. Their value is not checked.
    if constexpr (!ddc::in_tags_v<Patch, PatchSeqMin>) {
        EXPECT_NEAR(
                function_and_derivs(idx_dx, idx_xmin, idx_dy, idx_ymin),
                global_deriv_min_min,
                TOL);
        EXPECT_NEAR(
                function_and_derivs(idx_dx, idx_xmax, idx_dy, idx_ymin),
                global_deriv_max_min,
                TOL);
    }
    // For Patches in PatchSeqMax, we defined ddc::BoundCond::GREVILLE the local upper Y-boundary,
    // we don't need the cross-derivatives for y = ymax. Their value is not checked.
    if constexpr (!ddc::in_tags_v<Patch, PatchSeqMax>) {
        EXPECT_NEAR(
                function_and_derivs(idx_dx, idx_xmin, idx_dy, idx_ymax),
                global_deriv_min_max,
                TOL);
        EXPECT_NEAR(
                function_and_derivs(idx_dx, idx_xmax, idx_dy, idx_ymax),
                global_deriv_max_max,
                TOL);
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
        class... Patches,
        class... CoordTransformType>
void check_all_xy_derivatives(
        MultipatchField<DerivFieldOnPatch_host, Patches...>& functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        std::tuple<CoordTransformType...> coord_transforms,
        double const TOL = 5e-13)
{
    using PatchSeq = ddc::detail::TypeSeq<Patches...>;
    (check_xy_derivatives<Patches, PatchSeqMin, PatchSeqMax>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             function_g_coef,
             std::get<ddc::type_seq_rank_v<Patches, PatchSeq>>(coord_transforms),
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
        CoordTransformType const& coord_transform,
        double const TOL = 1e-14)
{
    using HostExecSpace = Kokkos::DefaultHostExecutionSpace;

    using DimX = typename Patch::Dim1;
    using DimY = typename Patch::Dim2;
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

    auto coord_transform_l_to_g = coord_transform.coord_transform.get_inverse_mapping();

    // Evaluate and compare the local and global spline representations ----------------------
    ddc::host_for_each(idx_range_xy, [&](typename Patch::Idx12 const idx) {
        Coord<Xg, Yg> eval_point_g(coord_transform_l_to_g(eval_points(idx)));
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
        class... Patches,
        class... CoordTransformType>
void check_all_spline_representation_agreement(
        MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        std::tuple<CoordTransformType...> coord_transforms,
        double const TOL = 5e-14)
{
    using PatchSeq = ddc::detail::TypeSeq<Patches...>;
    (check_spline_representation_agreement<Patches, PatchSeqMin, PatchSeqMax>(
             functions_and_derivs.template get<Patches>(),
             evaluator_g,
             get_const_field(function_g_coef),
             std::get<ddc::type_seq_rank_v<Patches, PatchSeq>>(coord_transforms),
             TOL),
     ...);
};
