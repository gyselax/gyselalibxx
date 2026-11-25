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


/**
 * @brief A class to evaluate all the interpolations of all the patches at once.
 *
 * This class computes the evaluations of the interpolations defined on each patch at
 * given coordinates, using the interpolation evaluators given for each patch.
 * The coordinates are stored in fields on each patch. On a given
 * storing patch, we do not enforce that a coordinate is physically located on this
 * storing patch. These fields of coordinates can be the result of characteristics
 * equations solving.
 *
 * This function is useful to avoid calling all the evaluators individually,
 * especially in a multipatch geometry with several patches. It also manages to
 * evaluate the right interpolation at a coordinate physically located on another patch
 * than the storing patch.
 *
 * Additionally, methods to compute the first derivatives and cross-derivatives
 * are implemented.
 *
 * Any 2D interpolation evaluator (see concepts::InterpolationEvaluator) can be used. Different
 * evaluator types can be used on different patches.
 *
 * @warning This operator does not work on batched domain.
 *
 * @warning This operator only works for patches defined on the same continuous dimensions!
 *
 * @tparam PatchSeq A ddc::TypeSeq of the patches.
 * @tparam EvaluatorSeq A ddc::TypeSeq of the 2D interpolation evaluators on each patch, in the same
 *                  order as the patches in PatchSeq.
 * @tparam PatchLocator A operator that finds the patch where a given coordinate is physically located.
 * @tparam ExtrapolationRule The extrapolation rule type for outside of the global domain.
 */
template <
        ddc::concepts::type_seq PatchSeq,
        ddc::concepts::type_seq EvaluatorSeq,
        class PatchLocator,
        class ExtrapolationRule>
class InterpolationEvaluator2DCollection;

/**
 * @brief A class to evaluate all the interpolations of all the patches at once.
 *
 * See the documentation of the primary template.
 *
 * @tparam Patches The patches.
 * @tparam Evaluators The 2D interpolation evaluators on each patch, in the same order as the patches.
 * @tparam PatchLocator A operator that finds the patch where a given coordinate is physically located.
 * @tparam ExtrapolationRule The extrapolation rule type for outside of the global domain.
 */
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

    /// A type alias to get the evaluator type on a specific patch.
    template <class Patch>
    using EvaluatorOnPatch = ddc::type_seq_element_t<
            ddc::type_seq_rank_v<Patch, PatchOrdering>,
            ddc::TypeSeq<Evaluators...>>;

    /// A type alias to get the traits of the evaluator on a specific patch.
    template <class Patch>
    using Traits = InterpolationEvaluatorTraits<EvaluatorOnPatch<Patch>>;

public:
    /// @brief The type of the 2D evaluation mesh used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using evaluation_idx_range_type = typename Traits<Patch>::evaluation_idx_range_type;

    /// @brief The type of the first discrete dimension of interest used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using evaluation_discrete_dimension_type1
            = ddc::type_seq_element_t<0, ddc::to_type_seq_t<evaluation_idx_range_type<Patch>>>;

    /// @brief The type of the second discrete dimension of interest used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using evaluation_discrete_dimension_type2
            = ddc::type_seq_element_t<1, ddc::to_type_seq_t<evaluation_idx_range_type<Patch>>>;

    /// @brief The type of the first evaluation continuous dimension used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using continuous_dimension_type1 =
            typename evaluation_discrete_dimension_type1<Patch>::continuous_dimension_type;

    /// @brief The type of the second evaluation continuous dimension used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using continuous_dimension_type2 =
            typename evaluation_discrete_dimension_type2<Patch>::continuous_dimension_type;

    /// @brief The discrete dimension representing the basis along first dimension.
    /// @tparam Patch Patch type.
    template <class Patch>
    using basis_type1 = ddc::
            type_seq_element_t<0, ddc::to_type_seq_t<typename Traits<Patch>::basis_idx_range_type>>;

    /// @brief The discrete dimension representing the basis along second dimension.
    /// @tparam Patch Patch type.
    template <class Patch>
    using basis_type2 = ddc::
            type_seq_element_t<1, ddc::to_type_seq_t<typename Traits<Patch>::basis_idx_range_type>>;


    /// @brief The type of the domain for the 1D evaluation mesh along first dimension used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using evaluation_idx_range_type1 = IdxRange<evaluation_discrete_dimension_type1<Patch>>;

    /// @brief The type of the domain for the 1D evaluation mesh along second dimension used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using evaluation_idx_range_type2 = IdxRange<evaluation_discrete_dimension_type2<Patch>>;

    /// @brief The type of the whole discrete domain (dimensions for grids and derivatives).
    /// @tparam ValuesField A Field or DerivField type storing the evaluated values.
    template <class ValuesField>
    using idx_range_type = typename ValuesField::discrete_domain_type;

    /// @brief TypeSeq of the derivative dimensions (contained in DerivField).
    /// @tparam ValuesField A Field or DerivField type storing the evaluated values.
    template <class ValuesField>
    using deriv_tags = detail::deriv_sub_set_t<ddc::to_type_seq_t<idx_range_type<ValuesField>>>;

    /// @brief TypeSeq of the grid dimensions (without derivative dimensions).
    /// @tparam ValuesField A Field or DerivField type storing the evaluated values.
    template <class ValuesField>
    using physical_grids = ddc::type_seq_remove_t<
            ddc::to_type_seq_t<idx_range_type<ValuesField>>,
            deriv_tags<ValuesField>>;

    /// @brief The type of the whole domain representing evaluation points.
    /// @tparam ValuesField A Field or DerivField type storing the evaluated values.
    template <class ValuesField>
    using batched_evaluation_idx_range_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<
                    physical_grids<ValuesField>>;

    /// @brief The type of the 2D index range on which the interpolation coefficients are defined.
    /// @tparam Patch Patch type.
    template <class Patch>
    using coeff_idx_range_type = typename Traits<Patch>::coeff_idx_range_type;

    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = typename FirstEvaluator::exec_space;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = typename FirstEvaluator::memory_space;

    /// @brief The data type that the data is saved on.
    using data_type = typename InterpolationEvaluatorTraits<FirstEvaluator>::data_type;

    /// @brief The ND coordinate type corresponding to the evaluation mesh.
    using coord_type = typename InterpolationEvaluatorTraits<FirstEvaluator>::coord_type;

    /**
     * @brief Type for MultipatchType: A field of 2D interpolation coefficients for a non-batched
     * interpolation defined on both of the Patch's logical dimensions.
     * Needed public for functions on GPU.
     */
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
    /// @brief The number of patches.
    static constexpr std::size_t n_patches = ddc::type_seq_size_v<PatchOrdering>;

    /// @brief The number of interpolation dimensions.
    static constexpr std::size_t rank()
    {
        return 2;
    }

private:
    // Asserts -----------------------------------------------------------------------------------
    /// Check that the values field can be used to store values evaluated on the StoringPatch.
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

    /// Check that the coefficients field can be used to evaluate the given values field.
    template <class StoringPatch, class ValuesField, class CoeffsField>
    static constexpr void check_coeffs_field_type()
    {
        static_assert(
                std::is_same_v<
                        typename CoeffsField::discrete_domain_type,
                        typename Traits<StoringPatch>::template batched_coeff_idx_range_type<
                                batched_evaluation_idx_range_type<ValuesField>>>,
                "The coefficients fields must be defined on the coefficient index range of the "
                "patch (with the same batch dimensions as the values field).");
        static_assert(
                std::is_same_v<typename CoeffsField::memory_space, memory_space>,
                "The coefficients fields need to be defined on the same memory space as the "
                "InterpolationEvaluator2DCollection class.");
        static_assert(
                std::is_const_v<typename CoeffsField::element_type>,
                "The coefficients should be constant.");
        static_assert(
                std::is_same_v<std::remove_const_t<typename CoeffsField::element_type>, data_type>,
                "The coefficients fields need to store the same data type as the "
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
    /**
     * @brief Instantiate a InterpolationEvaluator2DCollection.
     *
     * @warning The evaluators have to be sorted in the same order as the patches.
     *
     * @param evaluators The interpolation evaluators for each patch.
     * @param patch_locator An operator to locate a coordinate. The mapping stored in this class
     *      has to be invertible.
     * @param extrapolation_rule The extrapolation rules.
     */
    KOKKOS_FUNCTION InterpolationEvaluator2DCollection(
            Evaluators const&... evaluators,
            PatchLocator const& patch_locator,
            ExtrapolationRule const& extrapolation_rule)
        : m_patch_locator(patch_locator)
        , m_extrapolation_rule(extrapolation_rule)
        , m_evaluators(evaluators...)
    {
    }

    /// Default destructor callable from GPU
    KOKKOS_DEFAULTED_FUNCTION ~InterpolationEvaluator2DCollection() = default;

    /// Default copy operator callable from GPU
    KOKKOS_DEFAULTED_FUNCTION InterpolationEvaluator2DCollection(
            InterpolationEvaluator2DCollection const&)
            = default;

    // Assignment operators ----------------------------------------------------------------------
    /**
     * @brief Evaluate the 2D interpolations (described by their interpolation coefficients) at a given coordinate.
     *
     * The interpolation coefficients represent a 2D function defined on the basis of the evaluator of each patch.
     * They can be obtained via various methods, such as using an InterpolationBuilder2DCollection.
     * The coordinate is defined on one patch. But it can be physically located on another patch.
     *
     * @anchor InterpolationEvaluator2DCollectionOperator
     *
     * @param[in] coord_eval The coordinate where the interpolation is evaluated.
     * @param[in] patches_coeffs A MultipatchType of DField storing the 2D interpolation coefficients.
     * @return The value of the interpolation at the desired coordinate on the right patch.
     */
    template <class Coord>
    KOKKOS_FUNCTION data_type
    operator()(Coord const coord_eval, MultipatchCoeffs const& patches_coeffs) const
    {
        int const patch_idx = get_patch_idx(coord_eval);
        return recursive_dispatch_patch_function(Idx<>(), coord_eval, patches_coeffs, patch_idx);
    }

    /**
     * @brief Evaluate 2D interpolations (described by their interpolation coefficients) on meshes.
     *
     * See @ref InterpolationEvaluator2DCollectionOperator.
     *
     * @tparam MultipatchValues A MultipatchField of DField or DerivField.
     * @tparam MultipatchCoeffsType A MultipatchField of DField storing the 2D interpolation
     *          coefficients. The index range of the coefficients on each patch may contain
     *          batch dimensions (matching the batch dimensions of the values).
     * @param[out] patches_values A MultipatchType of DField or DerivField to store the values of the interpolations
     *          at the given coordinates.
     * @param[in] patches_coords A MultipatchType of Field of Coordinate storing the coordinates of the meshes.
     * @param[in] patches_coeffs A MultipatchType of DField storing the 2D interpolation coefficients.
     */
    template <class MultipatchValues, class MultipatchCoeffsType>
    void operator()(
            MultipatchValues const& patches_values,
            MultipatchCoordField const& patches_coords,
            MultipatchCoeffsType const& patches_coeffs) const
            requires(is_multipatch_field_v<MultipatchValues>)
    {
        static_assert(
                is_multipatch_field_v<MultipatchCoeffsType>,
                "The interpolation coefficients must be provided as a MultipatchField.");
        static_assert(
                std::is_same_v<typename MultipatchCoeffsType::PatchOrdering, PatchOrdering>,
                "The coefficients must be provided on the same patches as the values, "
                "in the same order.");
        (apply_evaluator_field<Patches>(
                 Idx<>(),
                 patches_values.template get<Patches>(),
                 patches_coords.template get<Patches>(),
                 patches_coeffs),
         ...);
    }


    // Derivatives operators ---------------------------------------------------------------------
    /**
     * @brief Differentiate 2D interpolations (described by their interpolation coefficients) at a given
     *  coordinate along a specified dimension of interest.
     *
     * See @ref InterpolationEvaluator2DCollectionOperator.
     * 
     * @warning The derivative cannot be computed outside of the domain. 
     *
     * @tparam DerivDims Derivative dimension types, it can be on the first dimension, the second
     * or both.
     * @param[in] deriv_order Order (k,l) of the derivatives (e.g. d_x^k d_y^l f).
     * @param[in] coord_eval The coordinate where the interpolation is differentiated.
     * @param[in] patches_coeffs A MultipatchType of DField storing the 2D interpolation coefficients.
     *
     * @return The derivative of the interpolation at the desired coordinate on the right patch.
     */
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

    /**
     * @brief Differentiate 2D interpolations (described by their interpolation coefficients) on meshes
     * along first dimension of interest.
     *
     * See @ref InterpolationEvaluator2DCollectionOperator.
     * 
     * @warning The derivatives cannot be computed outside of the domain. 
     *
     * @tparam DerivDims Derivative dimension types, it can be on the first dimension, the second
     * or both.
     * @param[in] deriv_order Order (k,l) of the derivatives (e.g. d_x^k d_y^l f(x,y)).
     * @param[out] patches_deriv A MultipatchType of DField or DerivField to store the derivatives of the interpolations
     *          at the given coordinates. 
     * @param[in] patches_coords A MultipatchType of Field of Coordinate storing the coordinates of the meshes.
     *          If patches_deriv is a DerivField, the given coordinates can be the coordinates for all function values. 
     *          The operator will select the correct coordinates. 
     * @param[in] patches_coeffs A MultipatchType of DField storing the 2D interpolation coefficients.
     *          The index range of the coefficients on each patch may contain batch dimensions
     *          (matching the batch dimensions of the values).
     * @tparam MultipatchValues A MultipatchField of DField or DerivField.
     * @tparam MultipatchCoeffsType A MultipatchField of DField storing the 2D interpolation coefficients.
     */
    template <class MultipatchValues, class MultipatchCoeffsType, class... DerivDims>
    void deriv(
            Idx<DerivDims...> const& deriv_order,
            MultipatchValues const& patches_deriv,
            MultipatchCoordField const& patches_coords,
            MultipatchCoeffsType const& patches_coeffs) const
            requires(is_multipatch_field_v<MultipatchValues>)
    {
        static_assert(
                (ddc::in_tags_v<DerivDims, ddc::TypeSeq<DerivDim1, DerivDim2>> && ...),
                "Please provide a deriv_order indexed on the derivative dimensions of the patches "
                "(same dimensions for all the patches).");
        static_assert(
                is_multipatch_field_v<MultipatchCoeffsType>,
                "The interpolation coefficients must be provided as a MultipatchField.");
        static_assert(
                std::is_same_v<typename MultipatchCoeffsType::PatchOrdering, PatchOrdering>,
                "The coefficients must be provided on the same patches as the values, "
                "in the same order.");
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
    /** @brief Integration of splines (described by their spline coefficients).
     *
     * See @ref InterpolationEvaluator2DCollectionOperator.
     *
     * This function is only available if all the evaluators are spline evaluators
     * (ddc::SplineEvaluatorND).
     *
     * @param[out] integrals An Kokkos::View on host containing the integrals of each spline on each patch.
     * @param[in] patches_coeffs A MultipatchType of DField storing the 2D interpolation coefficients.
     */
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

    /** @brief Compute the values or the derivatives of a given patch at the coordinates
     * defined on the given patch. Specialisation for patch_values on DField.
     * Needed public for functions on GPU. 
     * @tparam StoringPatch Patch type where the given coordinates are stored. 
     *      They are not especially physically located on this patch.
     * @tparam ValuesField The type of the field of values (DField or DerivField).
     * @tparam DerivDims Derivative dimension types, it can be on the first dimension, the second
     *      or both.
     * @param[in] deriv_order Order (k,l) of the derivatives (e.g. d_x^k d_y^l f).
     * @param[out] patch_values Field of values of the function or derivative.
     * @param[in] patch_coords ConstField of coordinates defined on the StoringPatch and
     *          where we want to evaluate the function or derivative.
     * @param[in] patches_coeffs MultipatchType of interpolation coefficients
     *          on every patches. The index range of the coefficients on each patch may
     *          contain batch dimensions (matching the batch dimensions of the values).
     */
    template <class StoringPatch, class ValuesField, class MultipatchCoeffsType, class... DerivDims>
    void apply_evaluator_field(
            Idx<DerivDims...> const& deriv_order,
            ValuesField const& patch_values,
            CoordConstFieldOnPatch<StoringPatch> const& patch_coords,
            MultipatchCoeffsType const patches_coeffs) const
    {
        check_values_field_type<StoringPatch, ValuesField>();
        check_coeffs_field_type<
                StoringPatch,
                ValuesField,
                typename MultipatchCoeffsType::template TypeOnPatch<StoringPatch>>();
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

    /** @brief Compute the derivatives of a given patch at the coordinates
     * defined on the given patch. Specialisation for patch_values on DerivField.
     * Needed public for functions on GPU. 
     * @tparam StoringPatch Patch type where the given coordinates are stored. 
     *      They are not especially physically located on this patch.
     * @tparam ValuesField The type of the DerivField of values.
     * @tparam DerivDims Derivative dimension types, it can be on the first dimension, the second
     * or both.
     * @param[in] deriv_order Order (k,l) of the derivatives (e.g. d_x^k d_y^l f).
     * @param[out] patch_values DerivField of values and derivatives of the function.
     * @param[in] patch_coords ConstField of coordinates defined on the StoringPatch and
     *          where we want to evaluate the function or derivative.
     * @param[in] patches_coeffs MultipatchType of interpolation coefficients
     *          on every patches. The index range of the coefficients on each patch may
     *          contain batch dimensions (matching the batch dimensions of the values).
     */
    template <class StoringPatch, class ValuesField, class MultipatchCoeffsType, class... DerivDims>
    void apply_evaluator_deriv_field(
            Idx<DerivDims...> const& deriv_order,
            ValuesField const& patch_values,
            CoordConstFieldOnPatch<StoringPatch> const& patch_coords,
            MultipatchCoeffsType const patches_coeffs) const
    {
        static_assert(is_deriv_field_v<ValuesField>, "The values field must be a DerivField.");
        check_values_field_type<StoringPatch, ValuesField>();
        check_coeffs_field_type<
                StoringPatch,
                ValuesField,
                typename MultipatchCoeffsType::template TypeOnPatch<StoringPatch>>();
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


    /** @brief Integrate the spline defined on the given patch.
     * @tparam Patch Patch type where the integration of the spline is computed. 
     * @param[out] integral Double, value of the integral of the spline on the given Patch.
     * @param[in] coeffs ConstField of spline coefficients of the spline on the given Patch.
     */
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

    /// @brief Dispatch the given coordinate on the right patch to evaluate the right interpolation.
    template <class IdxDeriv, class MultipatchCoeffsType, int TestPatchIdx = 0>
    KOKKOS_INLINE_FUNCTION data_type recursive_dispatch_patch_function(
            IdxDeriv const& deriv_order,
            Coord<CDim1, CDim2> coord,
            MultipatchCoeffsType const& patches_coeffs,
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
                        MultipatchCoeffsType,
                        TestPatchIdx + 1>(deriv_order, coord, patches_coeffs, patch_idx);
            }
        }
    }


    /**
     * @brief Call the evaluator of a patch with the given derivative order.
     *
     * The dimensions with a zero derivative order are removed from the derivative order
     * before calling the evaluator. If all the orders are zero the function itself is
     * evaluated. This is necessary as some evaluators (e.g. ddc::SplineEvaluatorND) do not
     * accept zero orders in deriv().
     */
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

    /// @brief Call the patch locator to get the index of the patch where the given coordinate
    /// is physically located.
    KOKKOS_INLINE_FUNCTION int get_patch_idx(Coord<CDim1, CDim2> const coord) const
    {
        using Mapping = typename PatchLocator::template get_mapping_on_logical_dim_t<CDim1, CDim2>;
        Mapping const mapping(m_patch_locator.template get_mapping_on_logical_dim<CDim1, CDim2>());
        return m_patch_locator(mapping(coord));
    }

    /// @brief Replace a coordinate inside the domain if it is periodic.
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
