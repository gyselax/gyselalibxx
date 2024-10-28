// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>
#include <tuple>
#include <utility>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"

/**
 * @brief A class to call all the builders of all the patches once.
 *
 * We need to instantiate all the builders for all the pacthes in the main code.
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
 * @tparam ExecSpace The space (CPU/GPU) where the calculations are carried out.
 * @tparam MemorySpace The space (CPU/GPU) where the coefficients and values are stored.
 * @tparam BSpline1OnPatch A type alias which provides the first BSpline type along which the splines are built.
 * @tparam BSpline2OnPatch A type alias which provides the second BSpline type along which the splines are built.
 * @tparam Grid1OnPatch A type alias which provides the first Grid type along which the interpolation points of the splines are found.
 * @tparam Grid2OnPatch A type alias which provides the second Grid type along which the interpolation points of the splines are found.
 * @tparam BcLower1 The lower boundary condition on the first dimension.
 * @tparam BcUpper1 The upper boundary condition on the first dimension.
 * @tparam BcLower2 The lower boundary condition on the second dimension.
 * @tparam BcUpper2 The upper boundary condition on the second dimension.
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
        ddc::BoundCond BcLower1,
        ddc::BoundCond BcUpper1,
        ddc::BoundCond BcLower2,
        ddc::BoundCond BcUpper2,
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
            "The BSpline1OnPatch argument does not define bsplines on the dimension where the "
            "grids defined by Grid1OnPatch are defined.");
    static_assert(
            ((std::is_same_v<
                    typename BSpline2OnPatch<Patches>::continuous_dimension_type,
                    typename Grid2OnPatch<Patches>::continuous_dimension_type>)&&...),
            "The BSpline2OnPatch argument does not define bsplines on the dimension where the "
            "grids defined by Grid2OnPatch are defined.");
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
        using type = ddc::SplineBuilder2D<
                ExecSpace,
                MemorySpace,
                BSpline1OnPatch<Patch>,
                BSpline2OnPatch<Patch>,
                Grid1OnPatch<Patch>,
                Grid2OnPatch<Patch>,
                BcLower1,
                BcUpper1,
                BcLower2,
                BcUpper2,
                Solver,
                Grid1D...>;
    };

    /// A type alias to get the builder type on a specific patch.
    template <class Patch>
    using BuilderOnPatch = typename Build_BuilderType<Patch, ValuesOnPatch<Patch>>::type;

    /// A type alias to get the batched spline coefficients on a specific patch.
    template <class Patch>
    using SplineOnPatch
            = DField<typename BuilderOnPatch<Patch>::batched_spline_domain_type, MemorySpace>;

    /// A type alias to get the batched derivatives along the first dimension on a specific patch.
    template <class Patch>
    using Derivs1OnPatch
            = DConstField<typename BuilderOnPatch<Patch>::batched_derivs_domain_type1, MemorySpace>;

    /// A type alias to get the batched derivatives along the first dimension on a specific patch.
    template <class Patch>
    using Derivs2OnPatch
            = DConstField<typename BuilderOnPatch<Patch>::batched_derivs_domain_type2, MemorySpace>;

    /// A type alias to get the batched cross-derivatives on a specific patch.
    template <class Patch>
    using Derivs12OnPatch
            = DConstField<typename BuilderOnPatch<Patch>::batched_derivs_domain_type, MemorySpace>;

    /// The type of the batched spline coefficients.
    using MultipatchSplineCoeffs = MultipatchField<SplineOnPatch, Patches...>;

    // For PERIODIC or GREVILLE boundary conditions
    /// The type of the values at the batched interpolation points.
    using MultipatchValues = MultipatchField<ValuesOnPatch, Patches...>;

    using MultipatchDerivs1 = MultipatchField<Derivs1OnPatch, Patches...>;

    using MultipatchDerivs2 = MultipatchField<Derivs2OnPatch, Patches...>;

    using MultipatchDerivs12 = MultipatchField<Derivs12OnPatch, Patches...>;

    /// The type of the internal storage of the SplineBuilders.
    using BuilderTuple = std::tuple<BuilderOnPatch<Patches> const&...>;


    BuilderTuple const m_builders;

private:
    template <class Patch, template <typename P> typename DerivTypeOnPatch>
    std::optional<DerivTypeOnPatch<Patch>> get_deriv_value(
            std::optional<MultipatchField<DerivTypeOnPatch, Patches...>> derivs) const
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
     * @param[in] mixed_derivs_min1_min2
     *      The values of the the cross-derivatives at the lower boundary in the first dimension
     *      and the lower boundary in the second dimension.
     * @param[in] mixed_derivs_max1_min2
     *      The values of the the cross-derivatives at the upper boundary in the first dimension
     *      and the lower boundary in the second dimension.
     * @param[in] mixed_derivs_min1_max2
     *      The values of the the cross-derivatives at the lower boundary in the first dimension
     *      and the upper boundary in the second dimension.
     * @param[in] mixed_derivs_max1_max2
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
            std::optional<MultipatchDerivs12> mixed_derivs_min1_min2 = std::nullopt,
            std::optional<MultipatchDerivs12> mixed_derivs_max1_min2 = std::nullopt,
            std::optional<MultipatchDerivs12> mixed_derivs_min1_max2 = std::nullopt,
            std::optional<MultipatchDerivs12> mixed_derivs_max1_max2 = std::nullopt) const
    {
        ((std::get<BuilderOnPatch<Patches> const&>(m_builders)(
                 splines.template get<Patches>(),
                 get_const_field(values.template get<Patches>()),
                 get_deriv_value<Patches, Derivs1OnPatch>(derivs_min1),
                 get_deriv_value<Patches, Derivs1OnPatch>(derivs_max1),
                 get_deriv_value<Patches, Derivs2OnPatch>(derivs_min2),
                 get_deriv_value<Patches, Derivs2OnPatch>(derivs_max2),
                 get_deriv_value<Patches, Derivs12OnPatch>(mixed_derivs_min1_min2),
                 get_deriv_value<Patches, Derivs12OnPatch>(mixed_derivs_max1_min2),
                 get_deriv_value<Patches, Derivs12OnPatch>(mixed_derivs_min1_max2),
                 get_deriv_value<Patches, Derivs12OnPatch>(mixed_derivs_max1_max2))),
         ...);
    };
};
