

# File interpolation\_builder\_collection.hpp

[**File List**](files.md) **>** [**interpolation**](dir_54175eb71bf3bdb1ca864f2a0d352472.md) **>** [**interpolation\_builder\_collection.hpp**](interpolation__builder__collection_8hpp.md)

[Go to the documentation of this file](interpolation__builder__collection_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <optional>
#include <tuple>
#include <type_traits>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"

template <ddc::concepts::type_seq PatchSeq, ddc::concepts::type_seq BuilderSeq>
class InterpolationBuilderCollection;

template <concepts::Patch... Patches, concepts::InterpolationBuilder1D... Builders>
class InterpolationBuilderCollection<ddc::TypeSeq<Patches...>, ddc::TypeSeq<Builders...>>
{
    static_assert(sizeof...(Patches) > 0, "At least one patch must be provided.");
    static_assert(
            sizeof...(Patches) == sizeof...(Builders),
            "There must be exactly one builder per patch.");
    static_assert(
            ((std::is_same_v<
                      interpolation_grid_t<Builders>,
                      typename Patches::
                              Grid1> || std::is_same_v<interpolation_grid_t<Builders>, typename Patches::Grid2>)&&...),
            "The builders must be provided in the same order as the patches and must "
            "interpolate on one of the grids (Grid1 or Grid2) of their patch.");

    using PatchOrdering = ddc::TypeSeq<Patches...>;

    using FirstBuilder = ddc::type_seq_element_t<0, ddc::TypeSeq<Builders...>>;

public:
    using exec_space = typename FirstBuilder::exec_space;

    using memory_space = typename FirstBuilder::memory_space;

    using data_type = typename InterpolationBuilderTraits<FirstBuilder>::data_type;

private:
    static_assert(
            (std::is_same_v<typename Builders::exec_space, exec_space> && ...),
            "All builders must use the same execution space.");
    static_assert(
            (std::is_same_v<typename Builders::memory_space, memory_space> && ...),
            "All builders must use the same memory space.");
    static_assert(
            (std::is_same_v<
                     typename InterpolationBuilderTraits<Builders>::data_type,
                     data_type> && ...),
            "All builders must use the same data type.");

    template <class Patch>
    using BuilderOnPatch = ddc::type_seq_element_t<
            ddc::type_seq_rank_v<Patch, PatchOrdering>,
            ddc::TypeSeq<Builders...>>;

    template <class Patch>
    using Traits = InterpolationBuilderTraits<BuilderOnPatch<Patch>>;

    template <class Patch, class ValuesType>
    using CoeffsOnPatch
            = Field<data_type,
                    typename Traits<Patch>::template batched_basis_idx_range_type<
                            typename ValuesType::discrete_domain_type>,
                    memory_space>;

    template <class Patch, class ValuesType>
    using DerivsOnPatch = ConstField<
            data_type,
            typename Traits<Patch>::template batched_derivs_idx_range_type<
                    typename ValuesType::discrete_domain_type>,
            memory_space>;

    template <class MultipatchValues>
    using MultipatchCoeffs = detail::MultipatchField<
            PatchOrdering,
            ddc::TypeSeq<CoeffsOnPatch<
                    Patches,
                    typename MultipatchValues::template TypeOnPatch<Patches>>...>>;

    template <class MultipatchValues>
    using MultipatchDerivs = detail::MultipatchField<
            PatchOrdering,
            ddc::TypeSeq<DerivsOnPatch<
                    Patches,
                    typename MultipatchValues::template TypeOnPatch<Patches>>...>>;

    template <class Patch>
    using IdxRangeInterpolationOnPatch = typename Traits<Patch>::interpolation_idx_range_type;

    template <class Patch>
    using IdxRangeCoeffOnPatch = typename Traits<Patch>::coeff_idx_range_type;

    template <class MultipatchIdxRangeBatched>
    using MultipatchIdxRangeBatchedBasis = detail::MultipatchType<
            PatchOrdering,
            ddc::TypeSeq<typename Traits<Patches>::template batched_basis_idx_range_type<
                    typename MultipatchIdxRangeBatched::template TypeOnPatch<Patches>>...>>;

    template <class MultipatchIdxRangeBatched>
    using MultipatchIdxRangeBatchedDeriv = detail::MultipatchType<
            PatchOrdering,
            ddc::TypeSeq<typename Traits<Patches>::template batched_derivs_idx_range_type<
                    typename MultipatchIdxRangeBatched::template TypeOnPatch<Patches>>...>>;

    using BuilderTuple = std::tuple<Builders const&...>;

    BuilderTuple const m_builders;

public:
    using interpolation_idx_range_type = MultipatchType<IdxRangeInterpolationOnPatch, Patches...>;

    using coeff_idx_range_type = MultipatchType<IdxRangeCoeffOnPatch, Patches...>;

    static constexpr std::size_t rank()
    {
        return 1;
    }

    template <class MultipatchIdxRangeBatched>
    using batched_basis_idx_range_type = MultipatchIdxRangeBatchedBasis<MultipatchIdxRangeBatched>;

    template <class MultipatchIdxRangeBatched>
    using batched_derivs_idx_range_type = MultipatchIdxRangeBatchedDeriv<MultipatchIdxRangeBatched>;

private:
    template <class Patch>
    BuilderOnPatch<Patch> const& get_builder() const
    {
        return std::get<ddc::type_seq_rank_v<Patch, PatchOrdering>>(m_builders);
    }

    template <class Patch, class MultipatchDerivsType>
    std::optional<typename MultipatchDerivsType::template TypeOnPatch<Patch>> get_deriv_value(
            std::optional<MultipatchDerivsType> derivs) const
    {
        if (derivs.has_value()) {
            return derivs->template get<Patch>();
        } else {
            return std::nullopt;
        }
    }

public:
    explicit InterpolationBuilderCollection(Builders const&... builders)
        : m_builders(std::tie(builders...)) {};


    ~InterpolationBuilderCollection() = default;

    template <class MultipatchValues>
    void operator()(
            MultipatchCoeffs<MultipatchValues> coeffs,
            MultipatchValues const& values,
            std::optional<MultipatchDerivs<MultipatchValues>> derivs_xmin = std::nullopt,
            std::optional<MultipatchDerivs<MultipatchValues>> derivs_xmax = std::nullopt) const
            requires(is_multipatch_field_v<MultipatchValues>)
    {
        (get_builder<Patches>()(
                 coeffs.template get<Patches>(),
                 get_const_field(values.template get<Patches>()),
                 get_deriv_value<Patches>(derivs_xmin),
                 get_deriv_value<Patches>(derivs_xmax)),
         ...);
    };

    template <class MultipatchIdxRangeBatched>
    MultipatchIdxRangeBatchedBasis<MultipatchIdxRangeBatched> const batched_basis_idx_ranges(
            MultipatchIdxRangeBatched idx_ranges) const
    {
        return MultipatchIdxRangeBatchedBasis<MultipatchIdxRangeBatched>(batched_basis_idx_range(
                get_builder<Patches>(),
                idx_ranges.template get<Patches>())...);
    };
};
```


