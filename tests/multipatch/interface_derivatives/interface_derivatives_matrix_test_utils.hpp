
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "coord_transformation_tools.hpp"
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

// -----------------------------------------------------------------------------------------------
// INITIALISATION OPERATORS ----------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------

/// @brief Initialise the function with f(x,y) = cos(2/3*pi*x)sin(y).
template <class GridXg, class GridYg, class Layout>
void initialise_2D_function(DField<IdxRange<GridXg, GridYg>, Kokkos::HostSpace, Layout> function)
{
    using Xg = typename GridXg::continuous_dimension_type;
    using Yg = typename GridYg::continuous_dimension_type;
    ddc::host_for_each(get_idx_range(function), [&](Idx<GridXg, GridYg> idx) {
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
    inverse_mapping_t<typename CoordTransformType::XYTransform> coord_transform_l_to_g_2d
            = coord_transform.coord_transform.get_inverse_mapping();
    ddc::host_for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        Coord<Xg, Yg> equiv_global_coord(coord_transform_l_to_g_2d(ddc::coordinate(idx)));
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

template <class Patch, class GridDeriv, class CoordTransformType, class SplineXYgEvaluator>
void initialise_derivatives(
        Idx<GridDeriv> idx_deriv_pos,
        DerivFieldOnPatch_host<Patch> function_and_derivs,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<
                typename SplineXYgEvaluator::bsplines_type1,
                typename SplineXYgEvaluator::bsplines_type2>>> const& const_function_g_coef,
        CoordTransformType coord_transform)
{
    using BSplinesXg = typename SplineXYgEvaluator::bsplines_type1;
    using BSplinesYg = typename SplineXYgEvaluator::bsplines_type2;
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    using DimDerivG = typename GridDeriv::continuous_dimension_type;
    static_assert(ddc::in_tags_v<DimDerivG, ddc::detail::TypeSeq<Xg, Yg>>);

    auto coord_transform_g_to_l_1d
            = get_1d_transform<std::is_same_v<DimDerivG, Xg>>(coord_transform);
    inverse_mapping_t<typename CoordTransformType::XYTransform> coord_transform_l_to_g_2d
            = coord_transform.coord_transform.get_inverse_mapping();

    using DerivCoord = decltype(coord_transform_g_to_l_1d)::CoordResult;
    using DimDeriv = ddc::type_seq_element_t<0, ddc::to_type_seq_t<DerivCoord>>;
    using GridAlongDeriv = std::conditional_t<
            std::is_same_v<typename Patch::Dim1, DimDeriv>,
            typename Patch::Grid1,
            typename Patch::Grid2>;
    using GridPerpToDeriv = std::conditional_t<
            std::is_same_v<typename Patch::Dim1, DimDeriv>,
            typename Patch::Grid2,
            typename Patch::Grid1>;

    Coord<DimDerivG> deriv_coord_glob = ddc::coordinate(idx_deriv_pos);
    Coord<DimDeriv> deriv_coord = coord_transform_g_to_l_1d(deriv_coord_glob);

    IdxRange<GridAlongDeriv> idx_range_derivs(get_idx_range(function_and_derivs));
    IdxRange<GridPerpToDeriv> idx_range_perp(get_idx_range(function_and_derivs));

    Idx<ddc::Deriv<DimDeriv>> idx_deriv(1);
    Idx<ddc::Deriv<DimDerivG>> idx_deriv_g(1);
    Idx<GridAlongDeriv> idx_deriv_pos_local
            = (std::fabs(deriv_coord - ddc::coordinate(idx_range_derivs.front())) < 1e-14)
                      ? idx_range_derivs.front()
                      : idx_range_derivs.back();

    ddc::host_for_each(idx_range_perp, [&](Idx<GridPerpToDeriv> const& idx_perp) {
        Coord<Xg, Yg>
                coord_glob(deriv_coord_glob, coord_transform_l_to_g_2d(ddc::coordinate(idx_perp)));
        function_and_derivs(idx_deriv, idx_deriv_pos_local, idx_perp)
                = coord_transform_g_to_l_1d.jacobian(deriv_coord_glob)
                  * evaluator_g.deriv(idx_deriv_g, coord_glob, const_function_g_coef);
    });
}

template <
        class Patch,
        class IdxRangeType,
        class CoordTransformType,
        class SplineXYgEvaluator,
        class BSplinesXg = typename SplineXYgEvaluator::bsplines_type1,
        class BSplinesYg = typename SplineXYgEvaluator::bsplines_type2>
std::tuple<double, double, double, double> get_cross_derivatives(
        IdxRangeType idx_range,
        SplineXYgEvaluator const& evaluator_g,
        host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const& function_g_coef,
        CoordTransformType const& coord_transform_handler)
{
    using Xg = typename BSplinesXg::continuous_dimension_type;
    using Yg = typename BSplinesYg::continuous_dimension_type;

    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;

    using IdxXY = Idx<GridX, GridY>;

    IdxRange<GridX> idx_range_x(idx_range);
    IdxRange<GridY> idx_range_y(idx_range);

    Idx<GridX> idx_xmin(idx_range_x.front());
    Idx<GridX> idx_xmax(idx_range_x.back());
    Idx<GridY> idx_ymin(idx_range_y.front());
    Idx<GridY> idx_ymax(idx_range_y.back());

    inverse_mapping_t<typename CoordTransformType::XYTransform> coord_transform_l_to_g_2d
            = coord_transform_handler.coord_transform.get_inverse_mapping();
    typename CoordTransformType::XTransform& coord_transform_g_to_l_1d_x
            = coord_transform_handler.x_transform;
    typename CoordTransformType::YTransform& coord_transform_g_to_l_1d_y
            = coord_transform_handler.y_transform;

    Coord<Xg, Yg> coord_min_min_g(
            coord_transform_l_to_g_2d(ddc::coordinate(IdxXY(idx_xmin, idx_ymin))));
    Coord<Xg, Yg> coord_max_min_g(
            coord_transform_l_to_g_2d(ddc::coordinate(IdxXY(idx_xmax, idx_ymin))));
    Coord<Xg, Yg> coord_min_max_g(
            coord_transform_l_to_g_2d(ddc::coordinate(IdxXY(idx_xmin, idx_ymax))));
    Coord<Xg, Yg> coord_max_max_g(
            coord_transform_l_to_g_2d(ddc::coordinate(IdxXY(idx_xmax, idx_ymax))));

    Idx<ddc::Deriv<Xg>, ddc::Deriv<Yg>> idx_deriv(1, 1);
    double global_deriv_min_min = coord_transform_g_to_l_1d_x.jacobian(Coord<Xg>(coord_min_min_g))
                                  * coord_transform_g_to_l_1d_y.jacobian(Coord<Yg>(coord_min_min_g))
                                  * evaluator_g.deriv(idx_deriv, coord_min_min_g, function_g_coef);
    double global_deriv_max_min = coord_transform_g_to_l_1d_x.jacobian(Coord<Xg>(coord_max_min_g))
                                  * coord_transform_g_to_l_1d_y.jacobian(Coord<Yg>(coord_max_min_g))
                                  * evaluator_g.deriv(idx_deriv, coord_max_min_g, function_g_coef);
    double global_deriv_min_max = coord_transform_g_to_l_1d_x.jacobian(Coord<Xg>(coord_min_max_g))
                                  * coord_transform_g_to_l_1d_y.jacobian(Coord<Yg>(coord_min_max_g))
                                  * evaluator_g.deriv(idx_deriv, coord_min_max_g, function_g_coef);
    double global_deriv_max_max = coord_transform_g_to_l_1d_x.jacobian(Coord<Xg>(coord_max_max_g))
                                  * coord_transform_g_to_l_1d_y.jacobian(Coord<Yg>(coord_max_max_g))
                                  * evaluator_g.deriv(idx_deriv, coord_max_max_g, function_g_coef);
    return std::make_tuple(
            global_deriv_min_min,
            global_deriv_max_min,
            global_deriv_min_max,
            global_deriv_max_max);
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
        CoordTransformType const& coord_transform_handler)
{
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;
    Idx<ddc::Deriv<typename Patch::Dim1>, ddc::Deriv<typename Patch::Dim2>> idx_dx_dy(1, 1);
    IdxRange<GridX, GridY> idx_range(get_idx_range(function_and_derivs));
    Idx<GridX> idx_xmin(idx_range.front());
    Idx<GridX> idx_xmax(idx_range.back());
    Idx<GridY> idx_ymin(idx_range.front());
    Idx<GridY> idx_ymax(idx_range.back());
    std::
            tie(function_and_derivs(idx_dx_dy, idx_xmin, idx_ymin),
                function_and_derivs(idx_dx_dy, idx_xmax, idx_ymin),
                function_and_derivs(idx_dx_dy, idx_xmin, idx_ymax),
                function_and_derivs(idx_dx_dy, idx_xmax, idx_ymax))
            = get_cross_derivatives<
                    Patch>(idx_range, evaluator_g, const_function_g_coef, coord_transform_handler);
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
        class DimDeriv,
        bool checkMin,
        bool checkMax,
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

    using DimDerivG = std::conditional_t<
            std::is_same_v<Coord<DimDeriv>, typename CoordTransformType::XTransform::CoordResult>,
            Xg,
            Yg>;

    using GridAlongDeriv = std::conditional_t<
            std::is_same_v<DimDeriv, typename Patch::Dim1>,
            typename Patch::Grid1,
            typename Patch::Grid2>;

    using GridPerpToDeriv = std::conditional_t<
            std::is_same_v<DimDeriv, typename Patch::Dim1>,
            typename Patch::Grid2,
            typename Patch::Grid1>;

    IdxRange<GridPerpToDeriv> idx_range_perp(get_idx_range(function_and_derivs));
    IdxRange<GridAlongDeriv> idx_range_par(get_idx_range(function_and_derivs));

    Idx<ddc::Deriv<DimDeriv>> idx_deriv(1);
    Idx<ddc::Deriv<DimDerivG>> idx_deriv_g(1);
    Idx<GridAlongDeriv> idx_par_min(idx_range_par.front());
    Idx<GridAlongDeriv> idx_par_max(idx_range_par.back());

    auto coord_transform_g_to_l_1d
            = get_1d_transform<std::is_same_v<DimDerivG, Xg>>(coord_transform);
    inverse_mapping_t<typename CoordTransformType::XYTransform> coord_transform_l_to_g_2d
            = coord_transform.coord_transform.get_inverse_mapping();

    ddc::host_for_each(idx_range_perp, [&](Idx<GridPerpToDeriv> const& idx_perp) {
        if constexpr (checkMin) {
            typename Patch::Idx12 idx_min(idx_par_min, idx_perp);
            Coord<Xg, Yg> interface_coord_min(coord_transform_l_to_g_2d(ddc::coordinate(idx_min)));
            double global_deriv_min
                    = coord_transform_g_to_l_1d.jacobian(Coord<DimDerivG>(interface_coord_min))
                      * evaluator_g.deriv(idx_deriv_g, interface_coord_min, function_g_coef);
            EXPECT_NEAR(
                    function_and_derivs(idx_deriv, idx_par_min, idx_perp),
                    global_deriv_min,
                    TOL);
        }

        if constexpr (checkMax) {
            typename Patch::Idx12 idx_max(idx_par_max, idx_perp);
            Coord<Xg, Yg> interface_coord_max(coord_transform_l_to_g_2d(ddc::coordinate(idx_max)));

            double global_deriv_max
                    = coord_transform_g_to_l_1d.jacobian(Coord<DimDerivG>(interface_coord_max))
                      * evaluator_g.deriv(idx_deriv_g, interface_coord_max, function_g_coef);
            EXPECT_NEAR(
                    function_and_derivs(idx_deriv, idx_par_max, idx_perp),
                    global_deriv_max,
                    TOL);
        }
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
    (check_derivatives<typename Patches::Dim1, true, true, Patches>(
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
    (check_derivatives<
             typename Patches::Dim2,
             !ddc::in_tags_v<Patches, PatchSeqMin>,
             !ddc::in_tags_v<Patches, PatchSeqMax>,
             Patches>(
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
    using GridX = typename Patch::Grid1;
    using GridY = typename Patch::Grid2;
    Idx<ddc::Deriv<typename Patch::Dim1>, ddc::Deriv<typename Patch::Dim2>> idx_dx_dy(1, 1);
    IdxRange<GridX, GridY> idx_range(get_idx_range(function_and_derivs));
    Idx<GridX> idx_xmin(idx_range.front());
    Idx<GridX> idx_xmax(idx_range.back());
    Idx<GridY> idx_ymin(idx_range.front());
    Idx<GridY> idx_ymax(idx_range.back());
    auto [global_deriv_min_min, global_deriv_max_min, global_deriv_min_max, global_deriv_max_max]
            = get_cross_derivatives<
                    Patch>(idx_range, evaluator_g, function_g_coef, coord_transform);

    // For Patches in PatchSeqMin, we defined ddc::BoundCond::GREVILLE the local lower Y-boundary,
    // we don't need the cross-derivatives for y = ymin. Their value is not checked.
    if constexpr (!ddc::in_tags_v<Patch, PatchSeqMin>) {
        EXPECT_NEAR(function_and_derivs(idx_dx_dy, idx_xmin, idx_ymin), global_deriv_min_min, TOL);
        EXPECT_NEAR(function_and_derivs(idx_dx_dy, idx_xmax, idx_ymin), global_deriv_max_min, TOL);
    }
    // For Patches in PatchSeqMax, we defined ddc::BoundCond::GREVILLE the local upper Y-boundary,
    // we don't need the cross-derivatives for y = ymax. Their value is not checked.
    if constexpr (!ddc::in_tags_v<Patch, PatchSeqMax>) {
        EXPECT_NEAR(function_and_derivs(idx_dx_dy, idx_xmin, idx_ymax), global_deriv_min_max, TOL);
        EXPECT_NEAR(function_and_derivs(idx_dx_dy, idx_xmax, idx_ymax), global_deriv_max_max, TOL);
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

    inverse_mapping_t<typename CoordTransformType::XYTransform> coord_transform_l_to_g_2d
            = coord_transform.coord_transform.get_inverse_mapping();

    // Evaluate and compare the local and global spline representations ----------------------
    ddc::host_for_each(idx_range_xy, [&](typename Patch::Idx12 const idx) {
        Coord<Xg, Yg> eval_point_g(coord_transform_l_to_g_2d(eval_points(idx)));
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
