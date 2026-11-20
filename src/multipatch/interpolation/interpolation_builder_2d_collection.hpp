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
#include "spline_builder_deriv_field_2d.hpp"
#include "types.hpp"

namespace detail {

/// A tool to check if a builder is a ddc::SplineBuilder2D.
template <class Builder>
inline constexpr bool is_ddc_spline_builder_2d_v = false;

/// Specialisation of is_ddc_spline_builder_2d_v for ddc::SplineBuilder2D.
template <
        class ExecSpace,
        class MemorySpace,
        class BSpline1,
        class BSpline2,
        class InterpolationDDim1,
        class InterpolationDDim2,
        ddc::SplineBuilderClosure BcLower1,
        ddc::SplineBuilderClosure BcUpper1,
        ddc::SplineBuilderClosure BcLower2,
        ddc::SplineBuilderClosure BcUpper2,
        ddc::SplineSolver Solver>
inline constexpr bool is_ddc_spline_builder_2d_v<ddc::SplineBuilder2D<
        ExecSpace,
        MemorySpace,
        BSpline1,
        BSpline2,
        InterpolationDDim1,
        InterpolationDDim2,
        BcLower1,
        BcUpper1,
        BcLower2,
        BcUpper2,
        Solver>> = true;

} // namespace detail

/**
 * @brief A class to call all the builders of all the patches once.
 *
 * We need to instantiate all the builders for all the patches in the main code.
 * We process the same way for the Field containing the interpolation coefficients and the
 * values of the function on each patch. The Fields are stored in MultipatchField
 * objects. The builders are stored in this class.
 * This class is instantiated with all the builders.
 * The operator() allows to call all the builders stored in the member of this class in
 * one single line.
 *
 * This function is useful to avoid calling all the builders individually, especially in
 * a multipatch geometry with several patches.
 *
 * Any 2D interpolation builder (see concepts::InterpolationBuilder) can be used. Different
 * builder types can be used on different patches (e.g. to use different boundary conditions).
 *
 * @warning InterpolationBuilder2DCollection is implemented for tensor-product multi-patch decompositions.
 *
 * @tparam PatchSeq A ddc::TypeSeq of the patches.
 * @tparam BuilderSeq A ddc::TypeSeq of the 2D interpolation builders on each patch, in the same
 *                  order as the patches in PatchSeq.
 */
template <ddc::concepts::type_seq PatchSeq, ddc::concepts::type_seq BuilderSeq>
class InterpolationBuilder2DCollection;

/**
 * @brief A class to call all the builders of all the patches once.
 *
 * See the documentation of the primary template.
 *
 * @tparam Patches The patches.
 * @tparam Builders The 2D interpolation builders on each patch, in the same order as the patches.
 */
template <concepts::Patch... Patches, concepts::InterpolationBuilder... Builders>
class InterpolationBuilder2DCollection<ddc::TypeSeq<Patches...>, ddc::TypeSeq<Builders...>>
{
    static_assert(sizeof...(Patches) > 0, "At least one patch must be provided.");
    static_assert(
            sizeof...(Patches) == sizeof...(Builders),
            "There must be exactly one builder per patch.");
    static_assert(
            ((InterpolationBuilderTraits<Builders>::rank() == 2) && ...),
            "All builders must be 2D builders.");

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

    static_assert(
            (std::is_same_v<
                     typename Traits<Patches>::interpolation_idx_range_type,
                     IdxRange<typename Patches::Grid1, typename Patches::Grid2>> && ...),
            "The builders must be provided in the same order as the patches and must "
            "interpolate on the grids (Grid1, Grid2) of their patch.");

    /// A structure to get the index range of the function values stored in a field.
    template <class FieldType>
    struct ValuesIdxRange
    {
        using type = typename FieldType::discrete_domain_type;
    };

    /// Specialisation of ValuesIdxRange for a field with derivatives.
    template <class FieldType>
    requires(is_deriv_field_v<FieldType>) struct ValuesIdxRange<FieldType>
    {
        using type = typename FieldType::physical_idx_range_type;
    };

    /// A type alias to get the batched interpolation coefficients on a specific patch.
    template <class Patch, class ValuesType>
    using CoeffsOnPatch
            = Field<data_type,
                    typename Traits<Patch>::template batched_basis_idx_range_type<
                            typename ValuesIdxRange<ValuesType>::type>,
                    memory_space>;

    /// A type alias to get the batched derivatives along the first dimension on a specific patch.
    template <class Patch, class ValuesType>
    using Derivs1OnPatch = ConstField<
            data_type,
            ddc::replace_dim_of_t<
                    typename ValuesType::discrete_domain_type,
                    typename Patch::Grid1,
                    ddc::Deriv<typename Patch::Dim1>>,
            memory_space>;

    /// A type alias to get the batched derivatives along the second dimension on a specific patch.
    template <class Patch, class ValuesType>
    using Derivs2OnPatch = ConstField<
            data_type,
            ddc::replace_dim_of_t<
                    typename ValuesType::discrete_domain_type,
                    typename Patch::Grid2,
                    ddc::Deriv<typename Patch::Dim2>>,
            memory_space>;

    /// A type alias to get the batched cross-derivatives index range on a specific patch.
    template <class Patch, class BatchedIdxRange>
    using batched_derivs_idx_range_type_on_patch = ddc::replace_dim_of_t<
            ddc::replace_dim_of_t<
                    BatchedIdxRange,
                    typename Patch::Grid1,
                    ddc::Deriv<typename Patch::Dim1>>,
            typename Patch::Grid2,
            ddc::Deriv<typename Patch::Dim2>>;

    /// A type alias to get the batched cross-derivatives on a specific patch.
    template <class Patch, class ValuesType>
    using Derivs12OnPatch = ConstField<
            data_type,
            batched_derivs_idx_range_type_on_patch<
                    Patch,
                    typename ValuesType::discrete_domain_type>,
            memory_space>;

    /// A type alias to get the DerivField type on a specific patch.
    template <class Patch>
    using DerivConstFieldOnPatch = DerivConstField<
            data_type,
            IdxRange<
                    ddc::Deriv<typename Patch::Dim1>,
                    typename Patch::Grid1,
                    ddc::Deriv<typename Patch::Dim2>,
                    typename Patch::Grid2>,
            memory_space>;

    /// A type alias to get the interpolation index range on a specific patch.
    template <class Patch>
    using IdxRangeInterpolationOnPatch = typename Traits<Patch>::interpolation_idx_range_type;

    /// A type alias to get the coefficient index range on a specific patch.
    template <class Patch>
    using IdxRangeCoeffOnPatch = typename Traits<Patch>::coeff_idx_range_type;

    /// The type of the batched interpolation coefficients.
    template <class MultipatchValues>
    using MultipatchCoeffs = detail::MultipatchField<
            PatchOrdering,
            ddc::TypeSeq<CoeffsOnPatch<
                    Patches,
                    typename MultipatchValues::template TypeOnPatch<Patches>>...>>;

    /// The type of the batched derivatives along the first dimension.
    template <class MultipatchValues>
    using MultipatchDerivs1 = detail::MultipatchField<
            PatchOrdering,
            ddc::TypeSeq<Derivs1OnPatch<
                    Patches,
                    typename MultipatchValues::template TypeOnPatch<Patches>>...>>;

    /// The type of the batched derivatives along the second dimension.
    template <class MultipatchValues>
    using MultipatchDerivs2 = detail::MultipatchField<
            PatchOrdering,
            ddc::TypeSeq<Derivs2OnPatch<
                    Patches,
                    typename MultipatchValues::template TypeOnPatch<Patches>>...>>;

    /// The type of the batched cross-derivatives.
    template <class MultipatchValues>
    using MultipatchDerivs12 = detail::MultipatchField<
            PatchOrdering,
            ddc::TypeSeq<Derivs12OnPatch<
                    Patches,
                    typename MultipatchValues::template TypeOnPatch<Patches>>...>>;

    /// The type of the function values and derivatives stored in DerivFields.
    using MultipatchDerivField = MultipatchField<DerivConstFieldOnPatch, Patches...>;

    template <class MultipatchIdxRangeBatched>
    using MultipatchIdxRangeBatchedBasis = detail::MultipatchType<
            PatchOrdering,
            ddc::TypeSeq<typename Traits<Patches>::template batched_basis_idx_range_type<
                    typename MultipatchIdxRangeBatched::template TypeOnPatch<Patches>>...>>;

    template <class MultipatchIdxRangeBatched>
    using MultipatchIdxRangeBatchedDeriv = detail::MultipatchType<
            PatchOrdering,
            ddc::TypeSeq<batched_derivs_idx_range_type_on_patch<
                    Patches,
                    typename MultipatchIdxRangeBatched::template TypeOnPatch<Patches>>...>>;

    /// The type of the internal storage of the builders.
    using BuilderTuple = std::tuple<Builders const&...>;

    BuilderTuple const m_builders;

public:
    /// @brief The ND index range for the interpolation mesh.
    using interpolation_idx_range_type = MultipatchType<IdxRangeInterpolationOnPatch, Patches...>;

    /// @brief The index range for the interpolation coefficients.
    using coeff_idx_range_type = MultipatchType<IdxRangeCoeffOnPatch, Patches...>;

    /// @brief The number of interpolation dimensions.
    static constexpr std::size_t rank()
    {
        return 2;
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
     * @brief Instantiate the InterpolationBuilder2DCollection from the builders on each patch.
     *
     * @warning The builders have to be sorted in the same order as the patches.
     *
     * @param builders Interpolation builders for each patch.
     */
    explicit InterpolationBuilder2DCollection(Builders const&... builders)
        : m_builders(std::tie(builders...)) {};


    ~InterpolationBuilder2DCollection() = default;

    /**
     * @brief Build the interpolation coefficients of each given function.
     *
     * @param[out] coeffs MultipatchField of all the Fields pointing to the interpolation coefficients.
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
     *
     * The derivatives are only used by builders which are ddc::SplineBuilder2D. They are ignored
     * by any other builder.
     *
     * @tparam MultipatchValues The type of the MultipatchField containing the function values.
     *          The index range of each field may contain batch dimensions.
     */
    template <class MultipatchValues>
    void operator()(
            MultipatchCoeffs<MultipatchValues> coeffs,
            MultipatchValues const& values,
            std::optional<MultipatchDerivs1<MultipatchValues>> derivs_min1 = std::nullopt,
            std::optional<MultipatchDerivs1<MultipatchValues>> derivs_max1 = std::nullopt,
            std::optional<MultipatchDerivs2<MultipatchValues>> derivs_min2 = std::nullopt,
            std::optional<MultipatchDerivs2<MultipatchValues>> derivs_max2 = std::nullopt,
            std::optional<MultipatchDerivs12<MultipatchValues>> cross_derivs_min1_min2
            = std::nullopt,
            std::optional<MultipatchDerivs12<MultipatchValues>> cross_derivs_max1_min2
            = std::nullopt,
            std::optional<MultipatchDerivs12<MultipatchValues>> cross_derivs_min1_max2
            = std::nullopt,
            std::optional<MultipatchDerivs12<MultipatchValues>> cross_derivs_max1_max2
            = std::nullopt) const
            requires(
                    is_multipatch_field_v<
                            MultipatchValues> && (!is_deriv_field_v<typename MultipatchValues::template TypeOnPatch<Patches>> && ...))
    {
        (apply_builder<Patches>(
                 coeffs.template get<Patches>(),
                 get_const_field(values.template get<Patches>()),
                 get_deriv_value<Patches>(derivs_min1),
                 get_deriv_value<Patches>(derivs_max1),
                 get_deriv_value<Patches>(derivs_min2),
                 get_deriv_value<Patches>(derivs_max2),
                 get_deriv_value<Patches>(cross_derivs_min1_min2),
                 get_deriv_value<Patches>(cross_derivs_max1_min2),
                 get_deriv_value<Patches>(cross_derivs_min1_max2),
                 get_deriv_value<Patches>(cross_derivs_max1_max2)),
         ...);
    };

    /**
     * @brief Build the interpolation coefficients of each given function.
     *
     * The derivatives are only used by builders which are ddc::SplineBuilder2D. Any other
     * builder only uses the function values.
     *
     * @param[out] coeffs MultipatchField of all the Fields pointing to the interpolation coefficients.
     * @param[in] functions_and_derivs MultipatchField of all the DerivFields pointing to the function values
     * and their derivatives.
     */
    void operator()(
            MultipatchCoeffs<MultipatchDerivField> coeffs,
            MultipatchDerivField functions_and_derivs) const
    {
        (apply_builder<Patches>(
                 coeffs.template get<Patches>(),
                 functions_and_derivs.template get<Patches>()),
         ...);
    };

    /**
     * @brief Get a MultipatchType collecting the index ranges typed on the whole basis dimensions.
     * These index ranges are adapted to allocate interpolation coefficients on the multi-patch domain.
     * @param[in] idx_ranges MultipatchType collecting the index ranges on the grids of the patches.
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

private:
    template <
            class Patch,
            class CoeffsField,
            class ValuesField,
            class Derivs1,
            class Derivs2,
            class Derivs12>
    void apply_builder(
            CoeffsField coeffs,
            ValuesField values,
            std::optional<Derivs1> derivs_min1,
            std::optional<Derivs1> derivs_max1,
            std::optional<Derivs2> derivs_min2,
            std::optional<Derivs2> derivs_max2,
            std::optional<Derivs12> cross_derivs_min1_min2,
            std::optional<Derivs12> cross_derivs_max1_min2,
            std::optional<Derivs12> cross_derivs_min1_max2,
            std::optional<Derivs12> cross_derivs_max1_max2) const
    {
        if constexpr (detail::is_ddc_spline_builder_2d_v<BuilderOnPatch<Patch>>) {
            get_builder<Patch>()(
                    coeffs,
                    values,
                    derivs_min1,
                    derivs_max1,
                    derivs_min2,
                    derivs_max2,
                    cross_derivs_min1_min2,
                    cross_derivs_max1_min2,
                    cross_derivs_min1_max2,
                    cross_derivs_max1_max2);
        } else {
            get_builder<Patch>()(coeffs, values);
        }
    };

    template <class Patch>
    void apply_builder(
            CoeffsOnPatch<Patch, DerivConstFieldOnPatch<Patch>> coeffs,
            DerivConstFieldOnPatch<Patch> function_and_deriv) const
    {
        using Builder = BuilderOnPatch<Patch>;
        if constexpr (detail::is_ddc_spline_builder_2d_v<Builder>) {
            using BasisSeq = ddc::to_type_seq_t<typename Traits<Patch>::coeff_idx_range_type>;
            SplineBuilderDerivField2D<
                    exec_space,
                    ddc::type_seq_element_t<0, BasisSeq>,
                    ddc::type_seq_element_t<1, BasisSeq>,
                    typename Patch::Grid1,
                    typename Patch::Grid2,
                    Builder::builder_type1::s_sbc_xmin,
                    Builder::builder_type1::s_sbc_xmax,
                    Builder::builder_type2::s_sbc_xmin,
                    Builder::builder_type2::s_sbc_xmax>
                    builder_applier(get_builder<Patch>());
            builder_applier(coeffs, function_and_deriv);
        } else {
            get_builder<Patch>()(coeffs, get_const_field(function_and_deriv.get_values_field()));
        }
    };
};
