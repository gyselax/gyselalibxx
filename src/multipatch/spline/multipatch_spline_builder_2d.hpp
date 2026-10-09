// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>
#include <tuple>
#include <utility>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"
#include "spline_builder_deriv_field_2d.hpp"

/**
 * @brief A class to call all the builders of all the patches once.
 *
 * We need to instantiate all the builders for all the patches in the main code.
 * We process the same way for the Field containing the spline coefficients and the
 * values of the function on each patch. The Fields are stored in MultipatchField
 * objects. The builders are stored in this class.
 * This class is instantiated with all the builders.
 * The operator() allows to call all the builders stored in the member of this class in 
 * one single line.
 *
 * This function is useful to avoid calling all the builders individually, especially in
 * a multipatch geometry with several patches.
 * 
 * @warning MultipatchSplineBuilder2D is implemented for tensor-product multi-patch decompositions. 
 *
 * @tparam ExecSpace The space (CPU/GPU) where the calculations are carried out.
 * @tparam MemorySpace The space (CPU/GPU) where the coefficients and values are stored.
 * @tparam BSpline1OnPatch A type alias which provides the first BSpline type along which the splines are built.
 * @tparam BSpline2OnPatch A type alias which provides the second BSpline type along which the splines are built.
 * @tparam Grid1OnPatch A type alias which provides the first Grid type along which the interpolation points of the splines are found.
 * @tparam Grid2OnPatch A type alias which provides the second Grid type along which the interpolation points of the splines are found.
 * @tparam SBCLower1 The lower spline closure on the first dimension.
 * @tparam SBCUpper1 The upper spline closure on the first dimension.
 * @tparam SBCLower2 The lower spline closure on the second dimension.
 * @tparam SBCUpper2 The upper spline closure on the second dimension.
 * @tparam BcTransition The boundary condition used at the interface between 2 patches.
 * @tparam Connectivity A MultipatchConnectivity object describing the interfaces between patches.
 * @tparam Solver The SplineSolver giving the backend used to perform the spline approximation. See DDC for more details.
 * @tparam ValuesOnPatch A type alias which provides the field type which will be used to pass the values of the function
 *                      at the interpolation points. The index range of this field type should contain any batch dimensions.
 */
template <
        class ExecSpace,
        class MemorySpace,
        template <typename P>
        typename BSpline1OnPatch,
        template <typename P>
        typename BSpline2OnPatch,
        template <typename P>
        typename Grid1OnPatch,
        template <typename P>
        typename Grid2OnPatch,
        ddc::SplineBuilderClosure SBCLower1,
        ddc::SplineBuilderClosure SBCUpper1,
        ddc::SplineBuilderClosure SBCLower2,
        ddc::SplineBuilderClosure SBCUpper2,
        ddc::SplineBuilderClosure BcTransition,
        class Connectivity,
        ddc::SplineSolver Solver,
        template <typename P>
        typename ValuesOnPatch,
        class... Patches>
class MultipatchSplineBuilder2D
{
    static_assert(
            (((ddc::is_uniform_bsplines_v<BSpline1OnPatch<Patches>>)
              || (ddc::is_non_uniform_bsplines_v<BSpline1OnPatch<Patches>>))
             && ...),
            "The BSpline1OnPatch argument does not create 1D BSpline objects.");
    static_assert(
            (((ddc::is_uniform_bsplines_v<BSpline2OnPatch<Patches>>)
              || (ddc::is_non_uniform_bsplines_v<BSpline2OnPatch<Patches>>))
             && ...),
            "The BSpline2OnPatch argument does not create 1D BSpline objects.");
    static_assert(
            (((ddc::is_uniform_point_sampling_v<Grid1OnPatch<Patches>>)
              || (ddc::is_non_uniform_point_sampling_v<Grid1OnPatch<Patches>>))
             && ...),
            "The Grid1OnPatch argument does not create 1D Grid objects.");
    static_assert(
            (((ddc::is_uniform_point_sampling_v<Grid2OnPatch<Patches>>)
              || (ddc::is_non_uniform_point_sampling_v<Grid2OnPatch<Patches>>))
             && ...),
            "The Grid2OnPatch argument does not create 1D Grid objects.");
    static_assert(
            ((std::is_same_v<
                    typename BSpline1OnPatch<Patches>::continuous_dimension_type,
                    typename Grid1OnPatch<Patches>::continuous_dimension_type>)&&...),
            "The BSpline1OnPatch argument does not define B-splines on the dimension where the "
            "grids defined by Grid1OnPatch are defined.");
    static_assert(
            ((std::is_same_v<
                    typename BSpline2OnPatch<Patches>::continuous_dimension_type,
                    typename Grid2OnPatch<Patches>::continuous_dimension_type>)&&...),
            "The BSpline2OnPatch argument does not define B-splines on the dimension where the "
            "grids defined by Grid2OnPatch are defined.");
    static_assert(
            !(is_deriv_field_v<ValuesOnPatch<Patches>> && ...),
            "ValuesOnPatch represents the field, not the field with derivatives.");

    /**
     * A small structure allowing the multiple grids to be unpacked from a field and repacked into
     * a SplineBuilder type.
     */
    template <class Patch, class FieldType>
    struct Build_BuilderType
    {
        static_assert(
                !std::is_same_v<Patch, Patch>,
                "The values should be saved in a constant field of doubles on the specified memory "
                "space.");
    };

    template <class Patch, class... Grid1D>
    struct Build_BuilderType<Patch, DConstField<IdxRange<Grid1D...>, MemorySpace>>
    {
        using lower_matching_edge1 = equivalent_edge_t<
                Edge<Patch, Grid1OnPatch<Patch>, FRONT>,
                typename Connectivity::interface_collection>;
        using upper_matching_edge1 = equivalent_edge_t<
                Edge<Patch, Grid1OnPatch<Patch>, BACK>,
                typename Connectivity::interface_collection>;
        using lower_matching_edge2 = equivalent_edge_t<
                Edge<Patch, Grid2OnPatch<Patch>, FRONT>,
                typename Connectivity::interface_collection>;
        using upper_matching_edge2 = equivalent_edge_t<
                Edge<Patch, Grid2OnPatch<Patch>, BACK>,
                typename Connectivity::interface_collection>;
        using type = ddc::SplineBuilder2D<
                ExecSpace,
                MemorySpace,
                BSpline1OnPatch<Patch>,
                BSpline2OnPatch<Patch>,
                Grid1OnPatch<Patch>,
                Grid2OnPatch<Patch>,
                std::is_same_v<lower_matching_edge1, OutsideEdge> ? SBCLower1 : BcTransition,
                std::is_same_v<upper_matching_edge1, OutsideEdge> ? SBCUpper1 : BcTransition,
                std::is_same_v<lower_matching_edge2, OutsideEdge> ? SBCLower2 : BcTransition,
                std::is_same_v<upper_matching_edge2, OutsideEdge> ? SBCUpper2 : BcTransition,
                Solver>;
    };

    /// A type alias to get the builder type on a specific patch.
    template <class Patch>
    using BuilderOnPatch = typename Build_BuilderType<Patch, ValuesOnPatch<Patch>>::type;

    /// A type alias to get the batched spline coefficients on a specific patch.
    template <class Patch>
    using SplineOnPatch = DField<
            typename BuilderOnPatch<Patch>::template batched_spline_domain_type<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    /// A type alias to get the batched derivatives along the first dimension on a specific patch.
    template <class Patch>
    using Derivs1OnPatch = DConstField<
            typename BuilderOnPatch<Patch>::template batched_derivs_domain_type1<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    /// A type alias to get the batched derivatives along the first dimension on a specific patch.
    template <class Patch>
    using Derivs2OnPatch = DConstField<
            typename BuilderOnPatch<Patch>::template batched_derivs_domain_type2<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    /// A type alias to get the batched cross-derivatives on a specific patch.
    template <class Patch>
    using Derivs12OnPatch = DConstField<
            typename BuilderOnPatch<Patch>::template batched_derivs_domain_type<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    /// A type alias to get the DerivField type on a specific patch.
    template <class Patch>
    using DerivConstFieldOnPatch = DerivField<
            const double,
            IdxRange<
                    ddc::Deriv<typename Patch::Dim1>,
                    typename Patch::Grid1,
                    ddc::Deriv<typename Patch::Dim2>,
                    typename Patch::Grid2>,
            MemorySpace>;

    /// A type alias to get the index range on the grids on a specific patch.
    template <class Patch>
    using IdxRangeOnPatch = typename Patch::IdxRange12;

    /// A type alias to get the index range on the B-splines on a specific patch.
    template <class Patch>
    using IdxRangeBSOnPatch = typename Patch::IdxRangeBS12;

    /// The type of the batched spline coefficients.
    using MultipatchSplineCoeffs = MultipatchField<SplineOnPatch, Patches...>;

    /// The type of the values at the batched interpolation points.
    using MultipatchValues = MultipatchField<ValuesOnPatch, Patches...>;

    using MultipatchDerivs1 = MultipatchField<Derivs1OnPatch, Patches...>;

    using MultipatchDerivs2 = MultipatchField<Derivs2OnPatch, Patches...>;

    using MultipatchDerivs12 = MultipatchField<Derivs12OnPatch, Patches...>;

    using MultipatchDerivField = MultipatchField<DerivConstFieldOnPatch, Patches...>;

    using MultipatchIdxRange = MultipatchType<IdxRangeOnPatch, Patches...>;

    using MultipatchIdxRangeBS = MultipatchType<IdxRangeBSOnPatch, Patches...>;

    /// The type of the internal storage of the SplineBuilders.
    using BuilderTuple = std::tuple<BuilderOnPatch<Patches> const&...>;


    BuilderTuple const m_builders;

private:
    template <class Patch, class MultipatchDerivs>
    std::optional<typename MultipatchDerivs::template TypeOnPatch<Patch>> get_deriv_value(
            std::optional<MultipatchDerivs> derivs) const
    {
        if (derivs.has_value()) {
            return derivs->template get<Patch>();
        } else {
            return std::nullopt;
        }
    }

public:
    /**
     * @brief Instantiate the MultipatchSplineBuilder from a std::tuple 
     * of all the builder on each patch. 
     * 
     * @warning The builders have to be sorted in the same order as the patches
     * in the tuple. 
     * 
     * @param builders Spline builders for each patch.
     */
    explicit MultipatchSplineBuilder2D(BuilderOnPatch<Patches> const&... builders)
        : m_builders(std::tie(builders...)) {};


    ~MultipatchSplineBuilder2D() = default;

    /**
     * @brief Build the spline representation of each given function.
     * 
     * @param[out] splines MultipatchField of all the Fields pointing to the spline representations. 
     * @param[in] values MultipatchField of all the Fields pointing to the function values. 
     * @param[in] derivs_min1 MultipatchField of all the ConstFields describing the function derivatives
     *                      in the first dimension at the lower bound of the second dimension.
     * @param[in] derivs_max1 MultipatchField of all the ConstFields describing the function derivatives
     *                      in the first dimension at the upper bound of the second dimension.
     * @param[in] derivs_min2 MultipatchField of all the ConstFields describing the function derivatives
     *                      in the second dimension at the lower bound of the first dimension.
     * @param[in] derivs_max2 MultipatchField of all the ConstFields describing the function derivatives
     *                      in the second dimension at the upper bound of the first dimension.
     * @param[in] cross_derivs_min1_min2
     *      The values of the the cross-derivatives at the lower boundary in the first dimension
     *      and the lower boundary in the second dimension.
     * @param[in] cross_derivs_max1_min2
     *      The values of the the cross-derivatives at the upper boundary in the first dimension
     *      and the lower boundary in the second dimension.
     * @param[in] cross_derivs_min1_max2
     *      The values of the the cross-derivatives at the lower boundary in the first dimension
     *      and the upper boundary in the second dimension.
     * @param[in] cross_derivs_max1_max2
     *      The values of the the cross-derivatives at the upper boundary in the first dimension
     *      and the upper boundary in the second dimension.
     */
    void operator()(
            MultipatchSplineCoeffs splines,
            MultipatchValues const& values,
            std::optional<MultipatchDerivs1> derivs_min1 = std::nullopt,
            std::optional<MultipatchDerivs1> derivs_max1 = std::nullopt,
            std::optional<MultipatchDerivs2> derivs_min2 = std::nullopt,
            std::optional<MultipatchDerivs2> derivs_max2 = std::nullopt,
            std::optional<MultipatchDerivs12> cross_derivs_min1_min2 = std::nullopt,
            std::optional<MultipatchDerivs12> cross_derivs_max1_min2 = std::nullopt,
            std::optional<MultipatchDerivs12> cross_derivs_min1_max2 = std::nullopt,
            std::optional<MultipatchDerivs12> cross_derivs_max1_max2 = std::nullopt) const
    {
        ((std::get<BuilderOnPatch<Patches> const&>(m_builders)(
                 splines.template get<Patches>(),
                 get_const_field(values.template get<Patches>()),
                 get_deriv_value<Patches>(derivs_min1),
                 get_deriv_value<Patches>(derivs_max1),
                 get_deriv_value<Patches>(derivs_min2),
                 get_deriv_value<Patches>(derivs_max2),
                 get_deriv_value<Patches>(cross_derivs_min1_min2),
                 get_deriv_value<Patches>(cross_derivs_max1_min2),
                 get_deriv_value<Patches>(cross_derivs_min1_max2),
                 get_deriv_value<Patches>(cross_derivs_max1_max2))),
         ...);
    };

    /**
     * @brief Build the spline representation of each given function.
     * 
     * @param[out] splines MultipatchField of all the Fields pointing to the spline representations. 
     * @param[in] functions_and_derivs MultipatchField of all the DerivFields pointing to the function values
     * and their derivatives. 
     */
    void operator()(MultipatchSplineCoeffs splines, MultipatchDerivField functions_and_derivs) const
    {
        (apply_builder<Patches>(
                 std::get<BuilderOnPatch<Patches> const&>(m_builders),
                 splines.template get<Patches>(),
                 functions_and_derivs.template get<Patches>()),
         ...);
    };

    /**
     * @brief Get a MultipatchType collecting the index ranges typed on the B-splines dimensions
     * of the patches. 
     * These index ranges are adapted to allocate splines on the multi-patch domain. 
     * @param[in] idx_ranges MultipatchType collecting the index ranges on the grids of the patches. 
     * @return a MultipatchType collecting the index ranges typed on the B-splines dimensions
     * of the patches. 
     */
    MultipatchIdxRangeBS const spline_idx_ranges(MultipatchIdxRange idx_ranges) const
    {
        return MultipatchIdxRangeBS(
                std::get<BuilderOnPatch<Patches> const&>(m_builders)
                        .batched_spline_domain(idx_ranges.template get<Patches>())...);
    };

private:
    template <class PatchP>
    static void apply_builder(
            BuilderOnPatch<PatchP> const& builder,
            SplineOnPatch<PatchP> spline,
            DerivConstFieldOnPatch<PatchP> function_and_deriv)
    {
        SplineBuilderDerivField2D<
                ExecSpace,
                BSpline1OnPatch<PatchP>,
                BSpline2OnPatch<PatchP>,
                Grid1OnPatch<PatchP>,
                Grid2OnPatch<PatchP>,
                BuilderOnPatch<PatchP>::builder_type1::s_sbc_xmin,
                BuilderOnPatch<PatchP>::builder_type1::s_sbc_xmax,
                BuilderOnPatch<PatchP>::builder_type2::s_sbc_xmin,
                BuilderOnPatch<PatchP>::builder_type2::s_sbc_xmax>
                builder_applier(builder);
        builder_applier(spline, function_and_deriv);
    };
};