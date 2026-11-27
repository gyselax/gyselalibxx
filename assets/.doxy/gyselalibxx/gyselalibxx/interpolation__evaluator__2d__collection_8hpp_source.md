

# File interpolation\_evaluator\_2d\_collection.hpp

[**File List**](files.md) **>** [**interpolation**](dir_54175eb71bf3bdb1ca864f2a0d352472.md) **>** [**interpolation\_evaluator\_2d\_collection.hpp**](interpolation__evaluator__2d__collection_8hpp.md)

[Go to the documentation of this file](interpolation__evaluator__2d__collection_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>
#include <utility>

#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "deriv_details.hpp"
#include "i_interpolation_evaluator.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "spline_interpolation.hpp"
#include "types.hpp"
#include "utils_patch_locators.hpp"
#include "view.hpp"


template <
        ddc::concepts::type_seq PatchSeq,
        ddc::concepts::type_seq EvaluatorSeq,
        class PatchLocator,
        class ExtrapolationRule>
class InterpolationEvaluator2DCollection;

template <
        concepts::Patch... Patches,
        concepts::InterpolationEvaluator... Evaluators,
        class PatchLocator,
        class ExtrapolationRule>
class InterpolationEvaluator2DCollection<
        ddc::TypeSeq<Patches...>,
        ddc::TypeSeq<Evaluators...>,
        PatchLocator,
        ExtrapolationRule>
{
    static_assert(sizeof...(Patches) > 0, "At least one patch must be provided.");
    static_assert(
            sizeof...(Patches) == sizeof...(Evaluators),
            "There must be exactly one evaluator per patch.");
    static_assert(
            ((InterpolationEvaluatorTraits<Evaluators>::rank() == 2) && ...),
            "All evaluators must be 2D evaluators.");

    // Patches
    using PatchOrdering = ddc::TypeSeq<Patches...>;

    using FirstEvaluator = ddc::type_seq_element_t<0, ddc::TypeSeq<Evaluators...>>;

public:
    template <class Patch>
    using EvaluatorOnPatch = ddc::type_seq_element_t<
            ddc::type_seq_rank_v<Patch, PatchOrdering>,
            ddc::TypeSeq<Evaluators...>>;

private:
    template <class Patch>
    using Traits = InterpolationEvaluatorTraits<EvaluatorOnPatch<Patch>>;

public:
    template <class Patch>
    using evaluation_idx_range_type = typename Traits<Patch>::evaluation_idx_range_type;

    template <class Patch>
    using evaluation_discrete_dimension_type1
            = ddc::type_seq_element_t<0, ddc::to_type_seq_t<evaluation_idx_range_type<Patch>>>;

    template <class Patch>
    using evaluation_discrete_dimension_type2
            = ddc::type_seq_element_t<1, ddc::to_type_seq_t<evaluation_idx_range_type<Patch>>>;

    template <class Patch>
    using continuous_dimension_type1 =
            typename evaluation_discrete_dimension_type1<Patch>::continuous_dimension_type;

    template <class Patch>
    using continuous_dimension_type2 =
            typename evaluation_discrete_dimension_type2<Patch>::continuous_dimension_type;

    template <class Patch>
    using basis_type1 = ddc::
            type_seq_element_t<0, ddc::to_type_seq_t<typename Traits<Patch>::basis_idx_range_type>>;

    template <class Patch>
    using basis_type2 = ddc::
            type_seq_element_t<1, ddc::to_type_seq_t<typename Traits<Patch>::basis_idx_range_type>>;


    template <class Patch>
    using evaluation_idx_range_type1 = IdxRange<evaluation_discrete_dimension_type1<Patch>>;

    template <class Patch>
    using evaluation_idx_range_type2 = IdxRange<evaluation_discrete_dimension_type2<Patch>>;

    template <class ValuesField>
    using idx_range_type = typename ValuesField::discrete_domain_type;

    template <class ValuesField>
    using deriv_tags = detail::deriv_sub_set_t<ddc::to_type_seq_t<idx_range_type<ValuesField>>>;

    template <class ValuesField>
    using physical_grids = ddc::type_seq_remove_t<
            ddc::to_type_seq_t<idx_range_type<ValuesField>>,
            deriv_tags<ValuesField>>;

    template <class ValuesField>
    using batched_evaluation_idx_range_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<
                    physical_grids<ValuesField>>;

    template <class Patch>
    using coeff_idx_range_type = typename Traits<Patch>::coeff_idx_range_type;

    using exec_space = typename FirstEvaluator::exec_space;

    using memory_space = typename FirstEvaluator::memory_space;

    using data_type = typename InterpolationEvaluatorTraits<FirstEvaluator>::data_type;

    using coord_type = typename InterpolationEvaluatorTraits<FirstEvaluator>::coord_type;

    template <class Patch>
    using CoeffOnPatch = ConstField<data_type, coeff_idx_range_type<Patch>, memory_space>;

private:
    static_assert(
            (std::is_same_v<typename Evaluators::exec_space, exec_space> && ...),
            "All evaluators must use the same execution space.");
    static_assert(
            (std::is_same_v<typename Evaluators::memory_space, memory_space> && ...),
            "All evaluators must use the same memory space.");
    static_assert(
            (std::is_same_v<
                     typename InterpolationEvaluatorTraits<Evaluators>::data_type,
                     data_type> && ...),
            "All evaluators must use the same data type.");
    static_assert(
            (std::is_same_v<
                     typename InterpolationEvaluatorTraits<Evaluators>::coord_type,
                     coord_type> && ...),
            "All evaluators must be defined on the same continuous dimensions.");
    static_assert(
            (std::is_same_v<
                     evaluation_idx_range_type<Patches>,
                     IdxRange<typename Patches::Grid1, typename Patches::Grid2>> && ...),
            "The evaluators must be provided in the same order as the patches and must "
            "evaluate on the grids (Grid1, Grid2) of their patch.");

    // Fields
    template <class Patch>
    using CoordConstFieldOnPatch
            = ConstField<coord_type, evaluation_idx_range_type<Patch>, memory_space>;

    // MultipatchTypes
    using MultipatchCoordField = MultipatchField<CoordConstFieldOnPatch, Patches...>;
    using MultipatchCoeffs = MultipatchField<CoeffOnPatch, Patches...>;
    using MultipatchEvaluators = detail::MultipatchType<PatchOrdering, ddc::TypeSeq<Evaluators...>>;


    // As we work with the same continuous dimensions --------------------------------------------
    // Continuous first dimension of the patches.
    using CDim1 = continuous_dimension_type1<ddc::type_seq_element_t<0, PatchOrdering>>;
    // Continuous second dimension of the patches.
    using CDim2 = continuous_dimension_type2<ddc::type_seq_element_t<0, PatchOrdering>>;

    // Dimension for the derivatives along the first dimension.
    using DerivDim1 = ddc::Deriv<CDim1>;
    // Dimension for the derivatives along the second dimension.
    using DerivDim2 = ddc::Deriv<CDim2>;

public:
    static constexpr std::size_t n_patches = ddc::type_seq_size_v<PatchOrdering>;

    static constexpr std::size_t rank()
    {
        return 2;
    }

private:
    // Asserts -----------------------------------------------------------------------------------
    template <class StoringPatch, class ValuesField>
    static constexpr void check_values_field_type()
    {
        static_assert(
                std::is_same_v<
                        evaluation_idx_range_type<StoringPatch>,
                        batched_evaluation_idx_range_type<ValuesField>>,
                "Operator not yet implemented for batched domains.");
        static_assert(
                std::is_same_v<typename ValuesField::memory_space, memory_space>,
                "The values fields need to be defined on the same memory space as the "
                "InterpolationEvaluator2DCollection class.");
    }

    static_assert(
            is_onion_patch_locator_v<PatchLocator>,
            "Operator currently works only on analytical mappings and patches on the same logical "
            "continuous dimension. E.g. OnionPatchLocator.");


    // Members -----------------------------------------------------------------------------------
    PatchLocator const m_patch_locator;
    ExtrapolationRule const m_extrapolation_rule;
    // The evaluators are stored by copy so they can be used on GPU.
    MultipatchEvaluators const m_evaluators;


public:
    KOKKOS_FUNCTION InterpolationEvaluator2DCollection(
            Evaluators const&... evaluators,
            PatchLocator const& patch_locator,
            ExtrapolationRule const& extrapolation_rule)
        : m_patch_locator(patch_locator)
        , m_extrapolation_rule(extrapolation_rule)
        , m_evaluators(evaluators...)
    {
    }

    KOKKOS_DEFAULTED_FUNCTION ~InterpolationEvaluator2DCollection() = default;

    KOKKOS_DEFAULTED_FUNCTION InterpolationEvaluator2DCollection(
            InterpolationEvaluator2DCollection const&)
            = default;

    template <concepts::Patch Patch>
    EvaluatorOnPatch<Patch> const& get() const
    {
        return m_evaluators.template get<Patch>();
    }

    // Assignment operators ----------------------------------------------------------------------
    template <class Coord>
    KOKKOS_FUNCTION data_type
    operator()(Coord const coord_eval, MultipatchCoeffs const& patches_coeffs) const
    {
        int const patch_idx = get_patch_idx(coord_eval);
        return recursive_dispatch_patch_function(Idx<>(), coord_eval, patches_coeffs, patch_idx);
    }

    template <class MultipatchValues>
    void operator()(
            MultipatchValues const& patches_values,
            MultipatchCoordField const& patches_coords,
            MultipatchCoeffs const& patches_coeffs) const
            requires(is_multipatch_field_v<MultipatchValues>)
    {
        (apply_evaluator_field<Patches>(
                 Idx<>(),
                 patches_values.template get<Patches>(),
                 patches_coords.template get<Patches>(),
                 patches_coeffs),
         ...);
    }


    // Derivatives operators ---------------------------------------------------------------------
    template <class... DerivDims>
    KOKKOS_FUNCTION data_type
    deriv(Idx<DerivDims...> const& deriv_order,
          Coord<CDim1, CDim2> const& coord_eval,
          MultipatchCoeffs const& patches_coeffs) const
    {
        static_assert(
                (ddc::in_tags_v<DerivDims, ddc::TypeSeq<DerivDim1, DerivDim2>> && ...),
                "Please provide a deriv_order indexed on the derivative dimensions of the patches "
                "(same dimensions for all the patches).");
        int const patch_idx = get_patch_idx(coord_eval);
        if (patch_idx < 0) {
            Kokkos::abort("The evaluation coordinate has to be on a patch."
                          "No extrapolation rule for derivatives. \n");
        }
        return recursive_dispatch_patch_function(
                deriv_order,
                coord_eval,
                patches_coeffs,
                patch_idx);
    }

    template <class MultipatchValues, class... DerivDims>
    void deriv(
            Idx<DerivDims...> const& deriv_order,
            MultipatchValues const& patches_deriv,
            MultipatchCoordField const& patches_coords,
            MultipatchCoeffs const& patches_coeffs) const
            requires(is_multipatch_field_v<MultipatchValues>)
    {
        static_assert(
                (ddc::in_tags_v<DerivDims, ddc::TypeSeq<DerivDim1, DerivDim2>> && ...),
                "Please provide a deriv_order indexed on the derivative dimensions of the patches "
                "(same dimensions for all the patches).");
        if constexpr ((is_deriv_field_v<
                               typename MultipatchValues::template TypeOnPatch<Patches>> && ...)) {
            (apply_evaluator_deriv_field<Patches>(
                     deriv_order,
                     patches_deriv.template get<Patches>(),
                     patches_coords.template get<Patches>(),
                     patches_coeffs),
             ...);
        } else {
            (apply_evaluator_field<Patches>(
                     deriv_order,
                     patches_deriv.template get<Patches>(),
                     patches_coords.template get<Patches>(),
                     patches_coeffs),
             ...);
        }
    }

    // Integrate operator ------------------------------------------------------------------------
    void integrate(
            DKokkosView_h<n_patches> const& integrals,
            MultipatchCoeffs const& patches_coeffs) const
            requires((detail::is_ddc_spline_evaluator_nd_v<Evaluators>)&&...)
    {
        (apply_integrate<Patches>(
                 integrals(ddc::type_seq_rank_v<Patches, PatchOrdering>),
                 patches_coeffs.template get<Patches>()),
         ...);
    }


public:
    // Apply functions to manage the values on each patch. ---------------------------------------

    template <class StoringPatch, class ValuesField, class... DerivDims>
    void apply_evaluator_field(
            Idx<DerivDims...> const& deriv_order,
            ValuesField const& patch_values,
            CoordConstFieldOnPatch<StoringPatch> const& patch_coords,
            MultipatchCoeffs const patches_coeffs) const
    {
        check_values_field_type<StoringPatch, ValuesField>();
        assert(get_idx_range(patch_values) == get_idx_range(patch_coords));

        using Index =
                typename batched_evaluation_idx_range_type<ValuesField>::discrete_element_type;
        batched_evaluation_idx_range_type<ValuesField> idx_range = get_idx_range(patch_values);

        const bool is_derivative = (deriv_order != order_zero<DerivDims...>());

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space(),
                idx_range,
                KOKKOS_CLASS_LAMBDA(Index const& idx) {
                    int const patch_idx = get_patch_idx(patch_coords(idx));
                    if (patch_idx < 0 && is_derivative) {
                        Kokkos::abort("The evaluation coordinate has to be on a patch."
                                      "No extrapolation rule for derivatives. \n");
                    }
                    patch_values(idx) = recursive_dispatch_patch_function(
                            deriv_order,
                            patch_coords(idx),
                            patches_coeffs,
                            patch_idx);
                });
    }

    template <class StoringPatch, class ValuesField, class... DerivDims>
    void apply_evaluator_deriv_field(
            Idx<DerivDims...> const& deriv_order,
            ValuesField const& patch_values,
            CoordConstFieldOnPatch<StoringPatch> const& patch_coords,
            MultipatchCoeffs const patches_coeffs) const
    {
        static_assert(is_deriv_field_v<ValuesField>, "The values field must be a DerivField.");
        check_values_field_type<StoringPatch, ValuesField>();
        if constexpr (ddc::in_tags_v<DerivDim1, ddc::TypeSeq<DerivDims...>>) {
            assert((IdxRange<DerivDim1>(patch_values.derivative_idx_range())
                            .contains(Idx<DerivDim1>(deriv_order)))
                   || (Idx<DerivDim1>(deriv_order) == Idx<DerivDim1>(0)));
        }
        if constexpr (ddc::in_tags_v<DerivDim2, ddc::TypeSeq<DerivDims...>>) {
            assert((IdxRange<DerivDim2>(patch_values.derivative_idx_range())
                            .contains(Idx<DerivDim2>(deriv_order)))
                   || (Idx<DerivDim2>(deriv_order) == Idx<DerivDim2>(0)));
        }

        using PGrid1 = typename StoringPatch::Grid1;
        using PGrid2 = typename StoringPatch::Grid2;

        IdxRangeSlice<PGrid1> idx_range_deriv_1
                = patch_values.template idx_range_for_deriv<PGrid1>();
        IdxRangeSlice<PGrid2> idx_range_deriv_2
                = patch_values.template idx_range_for_deriv<PGrid2>();

        IdxRange<PGrid1> idx_range_1(patch_values.idx_range());
        IdxRange<PGrid2> idx_range_2(patch_values.idx_range());

        // To match with the correct type.
        IdxRangeSlice<PGrid1> idx_range_to_slice_1(
                idx_range_1.front(),
                IdxStep<PGrid1>(idx_range_1.size()),
                IdxStep<PGrid1>(1));
        IdxRangeSlice<PGrid2> idx_range_to_slice_2(
                idx_range_2.front(),
                IdxStep<PGrid2>(idx_range_2.size()),
                IdxStep<PGrid2>(1));

        IdxRangeSlice<PGrid1, PGrid2> idx_range;
        if constexpr (std::is_same_v<Idx<DerivDims...>, Idx<DerivDim1>>) {
            idx_range = IdxRangeSlice<PGrid1, PGrid2>(idx_range_deriv_1, idx_range_to_slice_2);
        } else if constexpr (std::is_same_v<Idx<DerivDims...>, Idx<DerivDim2>>) {
            idx_range = IdxRangeSlice<PGrid1, PGrid2>(idx_range_to_slice_1, idx_range_deriv_2);
        } else if constexpr (std::is_same_v<Idx<DerivDims...>, Idx<DerivDim1, DerivDim2>>) {
            idx_range = IdxRangeSlice<PGrid1, PGrid2>(idx_range_deriv_1, idx_range_deriv_2);
        } else {
            idx_range = IdxRangeSlice<PGrid1, PGrid2>(idx_range_to_slice_1, idx_range_to_slice_2);
        }

        const bool is_derivative = (deriv_order != order_zero<DerivDims...>());

        // Correction if deriv_order is zero.
        if (!is_derivative) {
            idx_range = IdxRangeSlice<PGrid1, PGrid2>(idx_range_to_slice_1, idx_range_to_slice_2);
        }

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space(),
                idx_range,
                KOKKOS_CLASS_LAMBDA(Idx<PGrid1, PGrid2> const& idx) {
                    int const patch_idx = get_patch_idx(patch_coords(idx));
                    if (patch_idx < 0 && is_derivative) {
                        Kokkos::abort("The evaluation coordinate has to be on a patch."
                                      "No extrapolation rule for derivatives. \n");
                    }
                    patch_values(deriv_order, idx) = recursive_dispatch_patch_function(
                            deriv_order,
                            patch_coords(idx),
                            patches_coeffs,
                            patch_idx);
                });
    }


    template <class Patch>
    void apply_integrate(double& integral, CoeffOnPatch<Patch> const& coeffs) const
    {
        static_assert(
                std::is_same_v<exec_space, Kokkos::DefaultHostExecutionSpace>,
                "Can only be called on CPU: .integrals() not defined yet on GPU.");
        static_assert(
                Kokkos::SpaceAccessibility<exec_space, memory_space>::accessible,
                "Execution space and memory space have to be compatible.");

        using bsplines_1 = basis_type1<Patch>;
        using bsplines_2 = basis_type2<Patch>;
        using IdxBS12 = typename coeff_idx_range_type<Patch>::discrete_element_type;

        DFieldMem<IdxRange<bsplines_1>, memory_space> values1_alloc(
                "values1 (InterpolationEvaluator2DCollection::apply_integrate)",
                get_idx_range<bsplines_1>(coeffs));
        DField<IdxRange<bsplines_1>, memory_space> values1 = get_field(values1_alloc);
        DFieldMem<IdxRange<bsplines_2>, memory_space> values2_alloc(
                "values2 (InterpolationEvaluator2DCollection::apply_integrate)",
                get_idx_range<bsplines_2>(coeffs));
        DField<IdxRange<bsplines_2>, memory_space> values2 = get_field(values2_alloc);
        ddc::integrals(exec_space(), values1);
        ddc::integrals(exec_space(), values2);

        const std::source_location location = std::source_location::current();
        integral = ddc::parallel_transform_reduce(
                location.function_name(),
                exec_space(),
                get_idx_range(coeffs),
                0.0,
                ddc::reducer::sum<double>(),
                KOKKOS_LAMBDA(IdxBS12 const i12) {
                    return coeffs(i12) * values1(ddc::select<bsplines_1>(i12))
                           * values2(ddc::select<bsplines_2>(i12));
                });
    }



private:
    // Recursive method to dispatch the coordinates on the right patch ---------------------------

    template <class IdxDeriv, int TestPatchIdx = 0>
    KOKKOS_INLINE_FUNCTION data_type recursive_dispatch_patch_function(
            IdxDeriv const& deriv_order,
            Coord<CDim1, CDim2> coord,
            MultipatchCoeffs const& patches_coeffs,
            int const patch_idx) const
    {
        if constexpr (TestPatchIdx == ddc::type_seq_size_v<PatchOrdering>) {
            if (patch_idx >= 0) {
                Kokkos::abort("The recursion has reached the end without finding where "
                              "the coordinate is physically located.");
            }
            // Coord not on patch. Stop recursing.
            if constexpr (std::is_same_v<IdxDeriv, Idx<>>) {
                /* The operator currently works only for the case where the continuous 
                   dimensions of all the patches are the same. So the equivalent coordinates
                   are the same on any patches. 
                   Add a method to determine the patch from continuous dimensions, once 
                   we would like to deal with other cases. 
                */
                using AnyPatch = ddc::type_seq_element_t<0, PatchOrdering>;
                replace_periodic_coord_inside<AnyPatch>(coord);
                return m_extrapolation_rule(coord, patches_coeffs, patch_idx);
            } else {
                Kokkos::abort( // cppcheck-suppress missingReturn
                        "The derivatives cannot be evaluated at coordinates "
                        "outside of the domain.");
            }
        } else {
            if (patch_idx == TestPatchIdx) {
                using TestPatch = ddc::type_seq_element_t<TestPatchIdx, PatchOrdering>;
                replace_periodic_coord_inside<TestPatch>(coord);

                CoeffOnPatch<TestPatch> const test_coeffs
                        = patches_coeffs.template get<TestPatch>();
                return evaluate_on_patch(
                        m_evaluators.template get<TestPatch>(),
                        deriv_order,
                        coord,
                        test_coeffs);
            } else {
                return recursive_dispatch_patch_function<
                        IdxDeriv,
                        TestPatchIdx + 1>(deriv_order, coord, patches_coeffs, patch_idx);
            }
        }
    }


    template <class Evaluator, class Coeffs, class... DerivDims>
    KOKKOS_INLINE_FUNCTION data_type evaluate_on_patch(
            Evaluator const& evaluator,
            Idx<DerivDims...> const& deriv_order,
            Coord<CDim1, CDim2> const& coord,
            Coeffs const& coeffs) const
    {
        if constexpr (sizeof...(DerivDims) == 0) {
            return evaluator(coord, coeffs);
        } else if constexpr (sizeof...(DerivDims) == 1) {
            if (deriv_order == order_zero(deriv_order)) {
                return evaluator(coord, coeffs);
            } else {
                return evaluator.deriv(deriv_order, coord, coeffs);
            }
        } else {
            static_assert(sizeof...(DerivDims) == 2);
            Idx<DerivDim1> const order_1(deriv_order);
            Idx<DerivDim2> const order_2(deriv_order);
            bool const is_zero_1 = (order_1 == Idx<DerivDim1>(0));
            bool const is_zero_2 = (order_2 == Idx<DerivDim2>(0));
            if (is_zero_1 && is_zero_2) {
                return evaluator(coord, coeffs);
            } else if (is_zero_1) {
                return evaluator.deriv(order_2, coord, coeffs);
            } else if (is_zero_2) {
                return evaluator.deriv(order_1, coord, coeffs);
            } else {
                return evaluator.deriv(deriv_order, coord, coeffs);
            }
        }
    }

    // Useful functions --------------------------------------------------------------------------

    KOKKOS_INLINE_FUNCTION int get_patch_idx(Coord<CDim1, CDim2> const coord) const
    {
        using Mapping = typename PatchLocator::template get_mapping_on_logical_dim_t<CDim1, CDim2>;
        Mapping const mapping(m_patch_locator.template get_mapping_on_logical_dim<CDim1, CDim2>());
        return m_patch_locator(mapping(coord));
    }

    template <class Patch>
    KOKKOS_INLINE_FUNCTION void replace_periodic_coord_inside(Coord<CDim1, CDim2>& coord) const
    {
        using basis_1 = basis_type1<Patch>;
        using basis_2 = basis_type2<Patch>;

        Coord<CDim1> coord_1(coord);
        Coord<CDim2> coord_2(coord);

        ddcHelper::restrict_to_bspline_domain<basis_1>(coord_1);
        ddcHelper::restrict_to_bspline_domain<basis_2>(coord_2);

        coord = Coord<CDim1, CDim2>(coord_1, coord_2);
    }


    template <class... DerivDims>
    static KOKKOS_FUNCTION Idx<DerivDims...> order_zero()
    {
        return Idx<DerivDims...>(Idx<DerivDims>(0)...);
    }

    template <class... DerivDims>
    static KOKKOS_FUNCTION Idx<DerivDims...> order_zero(Idx<DerivDims...> const&)
    {
        return order_zero<DerivDims...>();
    }
};
```


