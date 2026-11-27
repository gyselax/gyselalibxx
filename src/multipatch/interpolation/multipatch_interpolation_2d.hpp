// SPDX-License-Identifier: MIT

#pragma once
#include <tuple>
#include <type_traits>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "i_interpolation.hpp"
#include "interpolation_builder_2d_collection.hpp"
#include "interpolation_evaluator_2d_collection.hpp"
#include "multipatch_type.hpp"
#include "type_seq_tools.hpp"
#include "types.hpp"

namespace detail {

/**
 * @brief An owning object that bundles the 2D interpolations of every patch of a
 * multipatch geometry.
 *
 * See the documentation of the alias MultipatchInterpolation2D, which should be used
 * instead of this class.
 *
 * This class can be used directly to define different interpolation types on each patch
 * (e.g. different boundary closures).
 *
 * @tparam PatchSeq A ddc::TypeSeq of the patches.
 * @tparam InterpolationSeq A ddc::TypeSeq of the 2D interpolations (see concepts::Interpolation)
 *          on each patch, in the order of PatchSeq.
 * @tparam PatchLocator The operator used to find the patch where a coordinate is physically located.
 * @tparam MultipatchExtrapolationRule The extrapolation rule used outside of the global domain.
 */
template <
        ddc::concepts::type_seq PatchSeq,
        ddc::concepts::type_seq InterpolationSeq,
        class PatchLocator,
        class MultipatchExtrapolationRule>
class MultipatchInterpolation2DImpl;

/**
 * @brief An owning object that bundles the 2D interpolations of every patch of a
 * multipatch geometry.
 *
 * See the documentation of the alias MultipatchInterpolation2D.
 *
 * @tparam Patches The patches.
 * @tparam Interpolations The 2D interpolations on each patch, in the same order as the patches.
 * @tparam PatchLocator The operator used to find the patch where a coordinate is physically located.
 * @tparam MultipatchExtrapolationRule The extrapolation rule used outside of the global domain.
 */
template <
        concepts::Patch... Patches,
        concepts::Interpolation... Interpolations,
        class PatchLocator,
        class MultipatchExtrapolationRule>
class MultipatchInterpolation2DImpl<
        ddc::TypeSeq<Patches...>,
        ddc::TypeSeq<Interpolations...>,
        PatchLocator,
        MultipatchExtrapolationRule>
{
    static_assert(sizeof...(Patches) > 0, "At least one patch must be provided.");
    static_assert(
            sizeof...(Patches) == sizeof...(Interpolations),
            "There must be exactly one interpolation per patch.");
    static_assert(
            ((Interpolations::rank() == 2) && ...),
            "The interpolations must be 2D interpolations on every patch.");

    using PatchOrdering = ddc::TypeSeq<Patches...>;

    template <class Patch>
    static constexpr std::size_t PatchIdx = ddc::type_seq_rank_v<Patch, PatchOrdering>;

    /// A type alias to get the interpolation type on a specific patch.
    template <class Patch>
    using InterpolationOnPatch
            = ddc::type_seq_element_t<PatchIdx<Patch>, ddc::TypeSeq<Interpolations...>>;

public:
    /// @brief The collection of the builders of each patch.
    using BuilderType = InterpolationBuilder2DCollection<
            PatchOrdering,
            ddc::TypeSeq<typename Interpolations::BuilderType...>>;

    /// @brief The collection of the evaluators of each patch.
    using EvaluatorType = InterpolationEvaluator2DCollection<
            PatchOrdering,
            ddc::TypeSeq<typename Interpolations::EvaluatorType...>,
            PatchLocator,
            MultipatchExtrapolationRule>;

    /// @brief The type of the index ranges on which the interpolations are built on each patch.
    using MultipatchIdxRange = ::MultipatchType<IdxRangeOnPatch, Patches...>;

    /// @brief The number of interpolation dimensions.
    static constexpr std::size_t rank()
    {
        return 2;
    }

private:
    // The interpolations must be declared before m_builder as m_builder stores
    // references to the builders owned by the interpolations.
    std::tuple<Interpolations...> m_interpolations;
    BuilderType m_builder;
    EvaluatorType m_evaluator;

public:
    /**
     * @brief Construct the interpolations on every patch.
     *
     * The interpolation of each patch is built from the index range of that patch.
     * The interpolations must therefore be constructible from an index range alone
     * (e.g. the extrapolation rules on each patch must be constructible automatically,
     * see is_extrapolation_rule_auto_constructible_v). If this is not the case, the
     * builders and evaluators should be built manually and collected with
     * InterpolationBuilder2DCollection and InterpolationEvaluator2DCollection.
     *
     * @param idx_ranges The 2D index ranges of the interpolation grids on each patch.
     * @param patch_locator An operator to locate a coordinate on the patches.
     * @param extrapolation_rule The extrapolation rule used outside of the global domain.
     */
    MultipatchInterpolation2DImpl(
            MultipatchIdxRange const& idx_ranges,
            PatchLocator const& patch_locator,
            MultipatchExtrapolationRule const& extrapolation_rule)
        : m_interpolations(idx_ranges.template get<Patches>()...)
        , m_builder(get_interpolation<Patches>().get_builder()...)
        , m_evaluator(
                  get_interpolation<Patches>().get_evaluator()...,
                  patch_locator,
                  extrapolation_rule)
    {
    }

    /**
     * @brief Construct the interpolations on every patch, using a default-constructed
     * extrapolation rule outside of the global domain.
     *
     * See the constructor which takes the extrapolation rule for more details.
     *
     * @param idx_ranges The 2D index ranges of the interpolation grids on each patch.
     * @param patch_locator An operator to locate a coordinate on the patches.
     */
    MultipatchInterpolation2DImpl(
            MultipatchIdxRange const& idx_ranges,
            PatchLocator const&
                    patch_locator) requires(std::
                                                    is_default_constructible_v<
                                                            MultipatchExtrapolationRule>)
        : MultipatchInterpolation2DImpl(idx_ranges, patch_locator, MultipatchExtrapolationRule())
    {
    }

    // The builder collection stores references to the builders in m_interpolations
    // so the object can be neither copied nor moved.
    MultipatchInterpolation2DImpl(MultipatchInterpolation2DImpl const&) = delete;
    MultipatchInterpolation2DImpl(MultipatchInterpolation2DImpl&&) = delete;
    MultipatchInterpolation2DImpl& operator=(MultipatchInterpolation2DImpl const&) = delete;
    MultipatchInterpolation2DImpl& operator=(MultipatchInterpolation2DImpl&&) = delete;

    ~MultipatchInterpolation2DImpl() = default;

    /**
     * @brief Return a const reference to the collection of the builders of each patch.
     * @return The BuilderType instance.
     */
    BuilderType const& get_builder() const
    {
        return m_builder;
    }

    /**
     * @brief Return a const reference to the collection of the evaluators of each patch.
     * @return The EvaluatorType instance.
     */
    EvaluatorType const& get_evaluator() const
    {
        return m_evaluator;
    }

    /**
     * @brief Return a const reference to the interpolation used on a given patch.
     * @tparam Patch The patch of interest.
     * @return The interpolation instance on the given patch.
     */
    template <class Patch>
    InterpolationOnPatch<Patch> const& get_interpolation() const
    {
        return std::get<PatchIdx<Patch>>(m_interpolations);
    }
};

/**
 * @brief The extrapolation rules used along a dimension on each patch by the aliases of
 * MultipatchInterpolation2D.
 *
 * The extrapolation on the patches is always superseded by the MultipatchExtrapolationRule,
 * so a null extrapolation is used. A periodic extrapolation is used along periodic dimensions,
 * as the interpolations require it.
 *
 * @tparam Dim The continuous dimension.
 */
template <class Dim>
using PatchExtrapolationRules = std::
        conditional_t<Dim::PERIODIC, ExtrapolationRule::Periodic, ExtrapolationRule::Null_Null>;

} // namespace detail

/**
 * @brief An owning object that bundles the 2D interpolations of every patch of a
 * multipatch geometry.
 *
 * MultipatchInterpolation2D constructs and owns one 2D interpolation (see
 * concepts::Interpolation, e.g. SplineInterpolator or LagrangeInterpolator) per patch,
 * and exposes them through an InterpolationBuilder2DCollection (get_builder()) and an
 * InterpolationEvaluator2DCollection (get_evaluator()). It is the multipatch
 * equivalent of the Interpolation classes and avoids building the builders, evaluators
 * and extrapolation rules of each patch manually.
 *
 * The patches are taken from the patch locator (PatchLocator::PatchOrdering).
 *
 * The builder collection refers to the builders owned by this object, so this object
 * must outlive any copy of the builder collection and can be neither copied nor moved.
 *
 * @tparam InterpOnPatch A type template such that InterpOnPatch<Patch> is the 2D
 *          interpolation used on the given patch. Its interpolation grids must be the
 *          grids (Grid1, Grid2) of the patch.
 * @tparam PatchLocator The operator used to find the patch where a coordinate is physically located.
 * @tparam MultipatchExtrapolationRule The extrapolation rule used outside of the global domain
 *          (e.g. ConstantExtrapolationRuleOnion or NullExtrapolationRule).
 */
template <
        template <class>
        class InterpOnPatch,
        class PatchLocator,
        class MultipatchExtrapolationRule>
using MultipatchInterpolation2D = detail::MultipatchInterpolation2DImpl<
        typename PatchLocator::PatchOrdering,
        type_seq_transform_t<InterpOnPatch, typename PatchLocator::PatchOrdering>,
        PatchLocator,
        MultipatchExtrapolationRule>;
