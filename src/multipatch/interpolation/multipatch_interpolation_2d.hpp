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

/// @brief The type of the builder of an interpolation.
template <class Interpolation>
using BuilderTypeOf = typename Interpolation::BuilderType;

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
 * @tparam MultipatchBuilder The multipatch builder. It must be constructible from a const
 *          reference to the InterpolationBuilder2DCollection of the builders of each patch,
 *          followed by any additional arguments passed to the constructor.
 */
template <
        ddc::concepts::type_seq PatchSeq,
        ddc::concepts::type_seq InterpolationSeq,
        class PatchLocator,
        class MultipatchExtrapolationRule,
        class MultipatchBuilder>
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
 * @tparam MultipatchBuilder The multipatch builder.
 */
template <
        concepts::Patch... Patches,
        concepts::Interpolation... Interpolations,
        class PatchLocator,
        class MultipatchExtrapolationRule,
        class MultipatchBuilder>
class MultipatchInterpolation2DImpl<
        ddc::TypeSeq<Patches...>,
        ddc::TypeSeq<Interpolations...>,
        PatchLocator,
        MultipatchExtrapolationRule,
        MultipatchBuilder>
{
    static_assert(sizeof...(Patches) > 0, "At least one patch must be provided.");
    static_assert(
            sizeof...(Patches) == sizeof...(Interpolations),
            "There must be exactly one interpolation per patch.");
    static_assert(
            ((Interpolations::rank() == 2) && ...),
            "The interpolations must be 2D interpolations on every patch.");

    using PatchOrdering = ddc::TypeSeq<Patches...>;

    using InterpolationSeq = ddc::TypeSeq<Interpolations...>;

    template <class Patch>
    static constexpr std::size_t PatchIdx = ddc::type_seq_rank_v<Patch, PatchOrdering>;

    /// A type alias to get the interpolation type on a specific patch.
    template <class Patch>
    using InterpolationOnPatch = ddc::type_seq_element_t<PatchIdx<Patch>, InterpolationSeq>;

public:
    /// @brief The collection of the (independent) builders of each patch.
    using LocalBuilderType = InterpolationBuilder2DCollection<
            PatchOrdering,
            type_seq_transform_t<BuilderTypeOf, InterpolationSeq>>;

    /// @brief The multipatch builder. This may be identical to LocalBuilderType but this
    /// is not an issue as InterpolationBuilder2DCollection only contains references.
    using BuilderType = MultipatchBuilder;

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
    // The order of the members matters: m_local_builder stores references to the builders
    // owned by m_interpolations and m_builder stores a reference to m_local_builder.
    // Similarly m_evaluator stores references to the evaluators owned by m_interpolations.
    std::tuple<Interpolations...> m_interpolations;
    LocalBuilderType m_local_builder;
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
     * @param builder_args Additional arguments passed to the constructor of the multipatch
     *          builder after the collection of the builders of each patch (e.g. an operator
     *          to estimate the derivatives at the patch interfaces).
     */
    template <class... BuilderArgs>
    MultipatchInterpolation2DImpl(
            MultipatchIdxRange const& idx_ranges,
            PatchLocator const& patch_locator,
            MultipatchExtrapolationRule const& extrapolation_rule,
            BuilderArgs&&... builder_args)
        : m_interpolations(idx_ranges.template get<Patches>()...)
        , m_local_builder(get_interpolation<Patches>().get_builder()...)
        , m_builder(m_local_builder, std::forward<BuilderArgs>(builder_args)...)
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

    // The builders store references to the builders in m_interpolations
    // so the object can be neither copied nor moved.
    MultipatchInterpolation2DImpl(MultipatchInterpolation2DImpl const&) = delete;
    MultipatchInterpolation2DImpl(MultipatchInterpolation2DImpl&&) = delete;
    MultipatchInterpolation2DImpl& operator=(MultipatchInterpolation2DImpl const&) = delete;
    MultipatchInterpolation2DImpl& operator=(MultipatchInterpolation2DImpl&&) = delete;

    ~MultipatchInterpolation2DImpl() = default;

    /**
     * @brief Return a const reference to the multipatch builder.
     * @return The BuilderType instance.
     */
    BuilderType const& get_builder() const
    {
        return m_builder;
    }

    /**
     * @brief Return a const reference to the collection of the (independent) builders of each patch.
     * @return The LocalBuilderType instance.
     */
    LocalBuilderType const& get_local_builder() const
    {
        return m_local_builder;
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
 * By default the builder (get_builder()) is the collection of the independent builders of
 * each patch (get_local_builder()). A different MultipatchBuilder (e.g. one which also
 * estimates the derivatives at the patch interfaces) can be provided. It is constructed
 * from the collection and any additional arguments passed to the constructor.
 *
 * The builders refer to the builders owned by this object, so this object must outlive
 * any copy of the builders and can be neither copied nor moved.
 *
 * @tparam InterpOnPatch A type template such that InterpOnPatch<Patch> is the 2D
 *          interpolation used on the given patch. Its interpolation grids must be the
 *          grids (Grid1, Grid2) of the patch.
 * @tparam PatchLocator The operator used to find the patch where a coordinate is physically located.
 * @tparam MultipatchExtrapolationRule The extrapolation rule used outside of the global domain
 *          (e.g. ConstantExtrapolationRuleOnion or NullExtrapolationRule).
 * @tparam MultipatchBuilder The multipatch builder (default: the InterpolationBuilder2DCollection
 *          of the builders of each patch).
 */
template <
        template <class>
        class InterpOnPatch,
        class PatchLocator,
        class MultipatchExtrapolationRule,
        class MultipatchBuilder = InterpolationBuilder2DCollection<
                typename PatchLocator::PatchOrdering,
                type_seq_transform_t<
                        detail::BuilderTypeOf,
                        type_seq_transform_t<InterpOnPatch, typename PatchLocator::PatchOrdering>>>>
using MultipatchInterpolation2D = detail::MultipatchInterpolation2DImpl<
        typename PatchLocator::PatchOrdering,
        type_seq_transform_t<InterpOnPatch, typename PatchLocator::PatchOrdering>,
        PatchLocator,
        MultipatchExtrapolationRule,
        MultipatchBuilder>;
