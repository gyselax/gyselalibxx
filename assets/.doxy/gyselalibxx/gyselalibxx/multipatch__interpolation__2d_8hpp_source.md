

# File multipatch\_interpolation\_2d.hpp

[**File List**](files.md) **>** [**interpolation**](dir_54175eb71bf3bdb1ca864f2a0d352472.md) **>** [**multipatch\_interpolation\_2d.hpp**](multipatch__interpolation__2d_8hpp.md)

[Go to the documentation of this file](multipatch__interpolation__2d_8hpp.md)


```C++
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

template <class Interpolation>
using BuilderTypeOf = typename Interpolation::BuilderType;

template <
        ddc::concepts::type_seq PatchSeq,
        ddc::concepts::type_seq InterpolationSeq,
        class PatchLocator,
        class MultipatchExtrapolationRule,
        class MultipatchBuilder>
class MultipatchInterpolation2DImpl;

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

    template <class Patch>
    using InterpolationOnPatch = ddc::type_seq_element_t<PatchIdx<Patch>, InterpolationSeq>;

public:
    using LocalBuilderType = InterpolationBuilder2DCollection<
            PatchOrdering,
            type_seq_transform_t<BuilderTypeOf, InterpolationSeq>>;

    using BuilderType = MultipatchBuilder;

    using EvaluatorType = InterpolationEvaluator2DCollection<
            PatchOrdering,
            ddc::TypeSeq<typename Interpolations::EvaluatorType...>,
            PatchLocator,
            MultipatchExtrapolationRule>;

    using MultipatchIdxRange = ::MultipatchType<IdxRangeOnPatch, Patches...>;

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

    BuilderType const& get_builder() const
    {
        return m_builder;
    }

    LocalBuilderType const& get_local_builder() const
    {
        return m_local_builder;
    }

    EvaluatorType const& get_evaluator() const
    {
        return m_evaluator;
    }

    template <class Patch>
    InterpolationOnPatch<Patch> const& get_interpolation() const
    {
        return std::get<PatchIdx<Patch>>(m_interpolations);
    }
};

template <class Dim>
using PatchExtrapolationRules = std::
        conditional_t<Dim::PERIODIC, ExtrapolationRule::Periodic, ExtrapolationRule::Null_Null>;

} // namespace detail

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
```


