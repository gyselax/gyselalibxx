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

/**
 * @brief A class to call all the builders of all the patches once.
 *
 * We need to instantiate all the builders for all the patches in the main code.
 * We process the same way for the Fields containing the interpolation coefficients and the
 * values of the function on each patch. The Fields are stored in MultipatchField
 * objects. The builders are stored in this class.
 * This class is instantiated with all the builders.
 * The operator() allows to call all the builders stored in the member of this class in
 * one single line.
 *
 * This function is useful to avoid calling all the builders individually, especially in
 * a multipatch geometry with several patches.
 *
 * Any 1D interpolation builder (see concepts::InterpolationBuilder1D) can be used. Different
 * builder types can be used on different patches (e.g. to use different boundary conditions).
 *
 * @tparam PatchSeq A ddc::TypeSeq of the patches.
 * @tparam BuilderSeq A ddc::TypeSeq of the 1D interpolation builders on each patch, in the same
 *                  order as the patches in PatchSeq.
 */
template <ddc::concepts::type_seq PatchSeq, ddc::concepts::type_seq BuilderSeq>
class InterpolationBuilderCollection;

/**
 * @brief A class to call all the builders of all the patches once.
 *
 * See the documentation of the primary template.
 *
 * @tparam Patches The patches.
 * @tparam Builders The 1D interpolation builders on each patch, in the same order as the patches.
 */
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
    /// @brief Execution space.
    using exec_space = typename FirstBuilder::exec_space;

    /// @brief Memory space.
    using memory_space = typename FirstBuilder::memory_space;

    /// @brief The data type that the data is saved on.
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

    /// A type alias to get the builder type on a specific patch.
    template <class Patch>
    using BuilderOnPatch = ddc::type_seq_element_t<
            ddc::type_seq_rank_v<Patch, PatchOrdering>,
            ddc::TypeSeq<Builders...>>;

    /// A type alias to get the traits of the builder on a specific patch.
    template <class Patch>
    using Traits = InterpolationBuilderTraits<BuilderOnPatch<Patch>>;

    /// A type alias to get the batched interpolation coefficients on a specific patch.
    template <class Patch, class ValuesType>
    using CoeffsOnPatch
            = Field<data_type,
                    typename Traits<Patch>::template batched_basis_idx_range_type<
                            typename ValuesType::discrete_domain_type>,
                    memory_space>;

    /// A type alias to get the batched derivatives on a specific patch.
    template <class Patch, class ValuesType>
    using DerivsOnPatch = ConstField<
            data_type,
            typename Traits<Patch>::template batched_derivs_idx_range_type<
                    typename ValuesType::discrete_domain_type>,
            memory_space>;

    /// The type of the batched interpolation coefficients.
    template <class MultipatchValues>
    using MultipatchCoeffs = detail::MultipatchField<
            PatchOrdering,
            ddc::TypeSeq<CoeffsOnPatch<
                    Patches,
                    typename MultipatchValues::template TypeOnPatch<Patches>>...>>;

    /// The type of the batched derivatives.
    template <class MultipatchValues>
    using MultipatchDerivs = detail::MultipatchField<
            PatchOrdering,
            ddc::TypeSeq<DerivsOnPatch<
                    Patches,
                    typename MultipatchValues::template TypeOnPatch<Patches>>...>>;

    /// A type alias to get the interpolation index range on a specific patch.
    template <class Patch>
    using IdxRangeInterpolationOnPatch = typename Traits<Patch>::interpolation_idx_range_type;

    /// A type alias to get the coefficient index range on a specific patch.
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

    /// The type of the internal storage of the builders.
    using BuilderTuple = std::tuple<Builders const&...>;

    BuilderTuple const m_builders;

public:
    /// @brief The 1D index range for the interpolation mesh.
    using interpolation_idx_range_type = MultipatchType<IdxRangeInterpolationOnPatch, Patches...>;

    /// @brief The index range for the interpolation coefficients.
    using coeff_idx_range_type = MultipatchType<IdxRangeCoeffOnPatch, Patches...>;

    /// @brief The number of interpolation dimensions.
    static constexpr std::size_t rank()
    {
        return 1;
    }

    /// @brief Multi-patch batched domain with the interpolation grids replaced by the basis grids.
    template <class MultipatchIdxRangeBatched>
    using batched_basis_idx_range_type = MultipatchIdxRangeBatchedBasis<MultipatchIdxRangeBatched>;

    /// @brief Multi-patch batched domain with the interpolation grids replaced by deriv types.
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
    /**
     * @brief Instantiate the InterpolationBuilderCollection from the builders on each patch.
     *
     * @warning The builders have to be sorted in the same order as the patches.
     *
     * @param builders Interpolation builders for each patch.
     */
    explicit InterpolationBuilderCollection(Builders const&... builders)
        : m_builders(std::tie(builders...)) {};


    ~InterpolationBuilderCollection() = default;

    /**
     * @brief Build the interpolation coefficients of each given function.
     *
     * @param coeffs MultipatchField of all the Fields pointing to the interpolation coefficients.
     * @param values MultipatchField of all the Fields pointing to the function values.
     * @param derivs_xmin MultipatchField of all the ConstFields describing the function derivatives at the lower bound.
     * @param derivs_xmax MultipatchField of all the ConstFields describing the function derivatives at the upper bound.
     *
     * @tparam MultipatchValues The type of the MultipatchField containing the function values.
     *          The index range of each field may contain batch dimensions.
     */
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

    /**
     * @brief Get a MultipatchType collecting the index ranges typed on the whole basis dimensions.
     * These index ranges are adapted to allocate interpolation coefficients on the multi-patch domain.
     * @param[in] idx_ranges MultipatchType collecting the (batched) index ranges on the grids of the patches.
     * @return a MultipatchType collecting the index ranges typed on the whole basis dimensions.
     */
    template <class MultipatchIdxRangeBatched>
    MultipatchIdxRangeBatchedBasis<MultipatchIdxRangeBatched> const batched_basis_idx_ranges(
            MultipatchIdxRangeBatched idx_ranges) const
    {
        return MultipatchIdxRangeBatchedBasis<MultipatchIdxRangeBatched>(batched_basis_idx_range(
                get_builder<Patches>(),
                idx_ranges.template get<Patches>())...);
    };
};
