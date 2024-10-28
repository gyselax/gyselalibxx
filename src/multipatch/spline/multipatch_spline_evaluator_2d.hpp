// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>
#include <utility>

#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "multipatch_type.hpp"
#include "types.hpp"
#include "utils_patch_locators.hpp"


/**
 * @brief A class to evaluate all the splines of all the patches at once.
 * 
 * This class computes the evaluations of the splines defined on each pacth at 
 * given coordinates. The class does not need to instantiate in advance individual 
 * spline evaluators. The coordinates are stored in fields on each patch. On a given 
 * storing patch, we do not enforce that a coordinate is physically located on this 
 * storing patch. These fields of coordinates can be the result of characteristics 
 * equations solving. 
 * 
 * This function is useful to avoid calling all the spline evaluators individually, 
 * especially in a multipatch geometry with several patches. It also manages to 
 * evaluate the right spline at a coordinate physically located on another patch 
 * than the storing patch.
 * 
 * Additionally, methods to compute the first derivatives and cross-derivatives 
 * are implemented. 
 * 
 * @warning This operator does not work on batched domain. 
 * 
 * @tparam ExecSpace The space (CPU/GPU) where the calculations are carried out.
 * @tparam MemorySpace The space (CPU/GPU) where the coefficients and values are stored.
 * @tparam BSpline1OnPatch A type alias which provides the first BSpline type along which the splines are built template on the Patch.
 * @tparam BSpline2OnPatch A type alias which provides the second BSpline type along which the splines are built template on the Patch.
 * @tparam Grid1OnPatch A type alias which provides the first Grid type along which the interpolation points 
 *          of the splines are found template on the Patch.
 * @tparam Grid2OnPatch A type alias which provides the second Grid type along which the interpolation points 
 *          of the splines are found template on the Patch.
 * @tparam ExtrapolationRule The extrapolation rule type for outside of the global domain. 
 * @tparam ValuesOnPatch A Field type storing the evaluated values of the splines. Template on the Patch. 
 * @tparam PatchLocator A operator that finds the patch where a given coordinate is physically located. 
 * @tparam Pacthes A variadic template of all the patches. 
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
        class ExtrapolationRule,
        template <typename P>
        typename ValuesOnPatch,
        class PatchLocator,
        class... Patches>
class MultipatchSplineEvaluator2D
{
public:
    /// @brief Tag to indicate that the value of the spline should be evaluated.
    struct eval_type
    {
    };

    /// @brief Tag to indicate that derivative of the spline should be evaluated.
    struct eval_deriv_type
    {
    };

    /// @brief The type of the first evaluation continuous dimension used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using continuous_dimension_type1 = typename BSpline1OnPatch<Patch>::continuous_dimension_type;

    /// @brief The type of the second evaluation continuous dimension used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using continuous_dimension_type2 = typename BSpline2OnPatch<Patch>::continuous_dimension_type;


    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = ExecSpace;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = MemorySpace;


    /// @brief The type of the first discrete dimension of interest used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using evaluation_discrete_dimension_type1 = Grid1OnPatch<Patch>;

    /// @brief The type of the second discrete dimension of interest used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using evaluation_discrete_dimension_type2 = Grid2OnPatch<Patch>;

    /// @brief The discrete dimension representing the B-splines along first dimension.
    /// @tparam Patch Patch type.
    template <class Patch>
    using bsplines_type1 = BSpline1OnPatch<Patch>;

    /// @brief The discrete dimension representing the B-splines along second dimension.
    /// @tparam Patch Patch type.
    template <class Patch>
    using bsplines_type2 = BSpline2OnPatch<Patch>;


    /// @brief The type of the domain for the 1D evaluation mesh along first dimension used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using evaluation_idx_range_type1 = IdxRange<evaluation_discrete_dimension_type1<Patch>>;

    /// @brief The type of the domain for the 1D evaluation mesh along second dimension used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using evaluation_idx_range_type2 = IdxRange<evaluation_discrete_dimension_type2<Patch>>;

    /// @brief The type of the domain for the 2D evaluation mesh used by this class.
    /// @tparam Patch Patch type.
    template <class Patch>
    using evaluation_idx_range_type = IdxRange<
            evaluation_discrete_dimension_type1<Patch>,
            evaluation_discrete_dimension_type2<Patch>>;


    /// @brief The type of the whole domain representing evaluation points.
    /// @tparam Patch Patch type.
    template <class Patch>
    using batched_evaluation_idx_range_type = typename ValuesOnPatch<Patch>::mdomain_type;


    /// @brief The type of the 1D spline domain corresponding to the first dimension of interest.
    /// @tparam Patch Patch type.
    template <class Patch>
    using spline_idx_range_type1 = IdxRange<bsplines_type1<Patch>>;

    /// @brief The type of the 1D spline domain corresponding to the second dimension of interest.
    /// @tparam Patch Patch type.
    template <class Patch>
    using spline_idx_range_type2 = IdxRange<bsplines_type2<Patch>>;

    /// @brief The type of the 2D spline domain corresponding to the dimensions of interest.
    /// @tparam Patch Patch type.
    template <class Patch>
    using spline_idx_range_type = IdxRange<bsplines_type1<Patch>, bsplines_type2<Patch>>;


private:
    template <class Patch>
    using CoordOnPatch
            = Coord<typename Grid1OnPatch<Patch>::continuous_dimension_type,
                    typename Grid2OnPatch<Patch>::continuous_dimension_type>;

    // Fields
    template <class Patch>
    using CoordConstFieldOnPatch
            = ConstField<CoordOnPatch<Patch>, evaluation_idx_range_type<Patch>, MemorySpace>;

public:
    /**
     * @brief Type for MultipatchType: A field of 2D spline coefficients for a non-batched spline defined
     * on both of the Patch's logical dimensions.
     * Needed public for functions on GPU. 
     */
    template <class Patch>
    using SplineCoeffOnPatch = DConstField<spline_idx_range_type<Patch>, MemorySpace>;

private:
    // MultipatchTypes
    using MultipatchValues = MultipatchField<ValuesOnPatch, Patches...>;
    using MultipatchCoordField = MultipatchField<CoordConstFieldOnPatch, Patches...>;
    using MultipatchSplineCoeff = MultipatchField<SplineCoeffOnPatch, Patches...>;

    // Patches
    using PatchOrdering = ddc::detail::TypeSeq<Patches...>;

public:
    /// @brief The number of patches.
    static constexpr std::size_t n_patches = ddc::type_seq_size_v<PatchOrdering>;

private:
    // Asserts -----------------------------------------------------------------------------------
    template <class Patch>
    static bool constexpr same_evaluation_idx_range_types = std::
            is_same_v<evaluation_idx_range_type<Patch>, batched_evaluation_idx_range_type<Patch>>;

    static_assert(
            (same_evaluation_idx_range_types<Patches> && ...),
            "Operator not yet implemented for batched domains.");

    static_assert(
            ((std::is_same_v<typename ValuesOnPatch<Patches>::memory_space, memory_space>)&&...),
            "Template ValuesOnPatch<Patch> need to be defined on the memory space than the"
            "MultipatchSplineEvaluator2D class.");

    static_assert(
            is_onion_patch_locator_v<PatchLocator>,
            "Operator currently works only on analytical mappings and patches on the same logical "
            "continous dimension. E.g. OnionPatchLocator.");


    // Members -----------------------------------------------------------------------------------
    PatchLocator const m_patch_locator;
    ExtrapolationRule const m_extrapolation_rule;


public:
    ~MultipatchSplineEvaluator2D() = default;

    /**
     * @brief Instantiate a MultipatchSplineEvaluator2D. 
     * 
     * @param patch_locator An operator to locate a coordinate. The mapping stored in this class 
     *      has to be invertible. 
     * @param extrapolation_rule The extrapolation rules.
     */
    MultipatchSplineEvaluator2D(
            PatchLocator const& patch_locator,
            ExtrapolationRule const& extrapolation_rule)
        : m_patch_locator(patch_locator)
        , m_extrapolation_rule(extrapolation_rule)
    {
    }


    // Assignment operators ----------------------------------------------------------------------
    /**
     * @brief Evaluate the 2D splines (described by their spline coefficients) at a given coordinate.
     *
     * The spline coefficients represent a 2D spline function defined on a B-splines (basis splines). 
     * They can be obtained via various methods, such as using a SplineBuilder2D or MultipatchSplineBuilder2D. 
     * The coordinate is defined on one patch. But it can be physically located on another patch. 
     * 
     * @anchor MultipatchSplineEvaluatorOperator
     * 
     * @tparam StoringPatch Patch type where the given coordinate is stored. It does not mean that 
     *      the coordinate is physically located on the patch. 
     * @param coord_eval The coordinate where the spline is evaluated. 
     * @param patches_splines A MultipatchType of DField storing the 2D spline coefficients.
     * @return The value of the spline at the desired coordinate on the right patch.
     */
    template <class Coord>
    KOKKOS_FUNCTION double operator()(
            Coord const coord_eval,
            MultipatchSplineCoeff const& patches_splines) const
    {
        int const patch_idx = get_patch_idx(coord_eval);
        return recursive_dispatch_patch_function<
                eval_type,
                eval_type>(coord_eval, patches_splines, patch_idx);
    }

    /**
     * @brief Evaluate 2D splines (described by their spline coefficients) on meshes.
     * 
     * See @ref MultipatchSplineEvaluatorOperator. 
     * 
     * @param[out] patches_values A MultipatchType of DField to store the values of the splines 
     *          at the given coordinates. 
     * @param[in] patches_coords A MultipatchType of Field of Coordinate storing the coordinates of the meshes.
     * @param[in] patches_splines A MultipatchType of DField storing the 2D spline coefficients.
     */
    void operator()(
            MultipatchValues const& patches_values,
            MultipatchCoordField const& patches_coords,
            MultipatchSplineCoeff const& patches_splines) const
    {
        (apply_evaluator<eval_type, eval_type, Patches>(
                 patches_values.template get<Patches>(),
                 patches_coords.template get<Patches>(),
                 patches_splines),
         ...);
    }


    // Derivatives operators ---------------------------------------------------------------------
    /**
     * @brief Differentiate 2D splines (described by their spline coefficients) at a given coordinate 
     * along first dimension of interest.
     *
     * See @ref MultipatchSplineEvaluatorOperator. 
     * 
     * @warning The derivative cannot be computed outside of the domain. 
     *      The coordinate do not still have to be defined on the right patch.
     *
     * @tparam StoringPatch Patch type where the given coordinate is stored. It does not mean that 
     *      the coordinate is physically located on the patch. 
     * @param coord_eval The coordinate where the spline is differentiated.
     * @param patches_splines A MultipatchType of DField storing the 2D spline coefficients.
     * 
     * @return The derivative of the spline at the desired coordinate on the right patch.
     */
    template <class Coord>
    KOKKOS_FUNCTION double deriv_dim_1(
            Coord const& coord_eval,
            MultipatchSplineCoeff const& patches_splines) const
    {
        int const patch_idx = get_patch_idx(coord_eval);
        if (patch_idx < 0) {
            Kokkos::abort("The evaluation coordinate has to be on a patch."
                          "No extrapolation rule for derivatives. \n");
        }
        return recursive_dispatch_patch_function<
                eval_deriv_type,
                eval_type>(coord_eval, patches_splines, patch_idx);
    }

    /**
     * @brief Differentiate 2D splines (described by their spline coefficients) at a given coordinate 
     * along second dimension of interest.
     *
     * See @ref MultipatchSplineEvaluatorOperator. 
     * 
     * @warning The derivative cannot be computed outside of the domain. 
     *
     * @tparam StoringPatch Patch type where the given coordinate is stored. It does not mean that 
     *      the coordinate is physically located on the patch. 
     * @param coord_eval The coordinate where the spline is differentiated.
     * @param patches_splines A MultipatchType of DField storing the 2D spline coefficients.
     * 
     * @return The derivative of the spline at the desired coordinate on the right patch.
     */
    template <class Coord>
    KOKKOS_FUNCTION double deriv_dim_2(
            Coord const& coord_eval,
            MultipatchSplineCoeff const& patches_splines) const
    {
        int const patch_idx = get_patch_idx(coord_eval);
        if (patch_idx < 0) {
            Kokkos::abort("The evaluation coordinate has to be on a patch."
                          "No extrapolation rule for derivatives. \n");
        }
        return recursive_dispatch_patch_function<
                eval_type,
                eval_deriv_type>(coord_eval, patches_splines, patch_idx);
    }

    /**
     * @brief Cross-differentiate 2D splines (described by their spline coefficients) at a given coordinate.
     *
     * See @ref MultipatchSplineEvaluatorOperator. 
     * 
     * @warning The derivative cannot be computed outside of the domain. 
     * 
     * @tparam StoringPatch Patch type where the given coordinate is stored. It does not mean that 
     *      the coordinate is physically located on the patch. 
     * @param coord_eval The coordinate where the spline is differentiated.
     * @param patches_splines A MultipatchType of DField storing the 2D spline coefficients.
     * 
     * @return The derivative of the spline at the desired coordinate on the right patch.
     */
    template <class Coord>
    KOKKOS_FUNCTION double deriv_1_and_2(
            Coord const& coord_eval,
            MultipatchSplineCoeff const& patches_splines) const
    {
        int const patch_idx = get_patch_idx(coord_eval);
        if (patch_idx < 0) {
            Kokkos::abort("The evaluation coordinate has to be on a patch."
                          "No extrapolation rule for derivatives. \n");
        }
        return recursive_dispatch_patch_function<
                eval_deriv_type,
                eval_deriv_type>(coord_eval, patches_splines, patch_idx);
        ;
    }


    /**
     * @brief Differentiate 2D splines (described by their spline coefficients) at a given coordinate 
     *  along a specified dimension of interest.
     *
     * See @ref MultipatchSplineEvaluatorOperator. 
     * 
     * @warning The derivative cannot be computed outside of the domain. 
     *
     * @tparam InterestDim Dimension of StoringPatch along which differentiation is performed.
     * @tparam StoringPatch Patch type where the given coordinate is stored. It does not mean that 
     *      the coordinate is physically located on the patch. 
     * @param coord_eval The coordinate where the spline is differentiated.
     * @param patches_splines A MultipatchType of DField storing the 2D spline coefficients.
     * 
     * @return The derivative of the spline at the desired coordinate on the right patch.
     */
    template <class InterestDim, class Dim1, class Dim2>
    KOKKOS_FUNCTION double deriv(
            Coord<Dim1, Dim2> const& coord_eval,
            MultipatchSplineCoeff const& patches_splines) const
    {
        static_assert((std::is_same_v<InterestDim, Dim1>) || (std::is_same_v<InterestDim, Dim2>));
        if constexpr (std::is_same_v<InterestDim, Dim1>) {
            return deriv_dim_1(coord_eval, patches_splines);
        } else if constexpr (std::is_same_v<InterestDim, Dim2>) {
            return deriv_dim_2(coord_eval, patches_splines);
        }
    }

    /**
     * @brief Differentiate 2D splines (described by their spline coefficients) on a meshes
     * along first dimension of interest.
     *
     * See @ref MultipatchSplineEvaluatorOperator. 
     * 
     * @warning The derivatives cannot be computed outside of the domain. 
     *
     * @param[out] patches_deriv_1 A MultipatchType of DField to store the derivatives of the splines 
     *          at the given coordinates. 
     * @param[in] patches_coords A MultipatchType of Field of Coordinate storing the coordinates of the meshes.
     * @param[in] patches_splines A MultipatchType of DField storing the 2D spline coefficients.
     */
    void deriv_dim_1(
            MultipatchValues const& patches_deriv_1,
            MultipatchCoordField const& patches_coords,
            MultipatchSplineCoeff const& patches_splines) const
    {
        (apply_evaluator<eval_deriv_type, eval_type, Patches>(
                 patches_deriv_1.template get<Patches>(),
                 patches_coords.template get<Patches>(),
                 patches_splines),
         ...);
    }

    /**
     * @brief Differentiate 2D splines (described by their spline coefficients) on a meshes
     * along second dimension of interest.
     *
     * See @ref MultipatchSplineEvaluatorOperator. 
     * 
     * @warning The derivatives cannot be computed outside of the domain. 
     *
     * @param[out] patches_deriv_2 A MultipatchType of DField to store the derivatives of the splines 
     *          at the given coordinates. 
     * @param[in] patches_coords A MultipatchType of Field of Coordinate storing the coordinates of the meshes.
     * @param[in] patches_splines A MultipatchType of DField storing the 2D spline coefficients.
     */
    void deriv_dim_2(
            MultipatchValues const& patches_deriv_2,
            MultipatchCoordField const& patches_coords,
            MultipatchSplineCoeff const& patches_splines) const
    {
        (apply_evaluator<eval_type, eval_deriv_type, Patches>(
                 patches_deriv_2.template get<Patches>(),
                 patches_coords.template get<Patches>(),
                 patches_splines),
         ...);
    }

    /** @brief Cross-differentiate 2D splines (described by their spline coefficients) on a meshes.
     * 
     * See @ref MultipatchSplineEvaluatorOperator. 
     * 
     * @warning The derivatives cannot be computed outside of the domain. 
     *
     * @param[out] patches_deriv_12 A MultipatchType of DField to store the cross-derivatives of the splines 
     *          at the given coordinates. 
     * @param[in] patches_coords A MultipatchType of Field of Coordinate storing the coordinates of the meshes.
     * @param[in] patches_splines A MultipatchType of DField storing the 2D spline coefficients.
     */
    void deriv_1_and_2(
            MultipatchValues const& patches_deriv_12,
            MultipatchCoordField const& patches_coords,
            MultipatchSplineCoeff const& patches_splines) const
    {
        (apply_evaluator<eval_deriv_type, eval_deriv_type, Patches>(
                 patches_deriv_12.template get<Patches>(),
                 patches_coords.template get<Patches>(),
                 patches_splines),
         ...);
    }


    // Integrate operator ------------------------------------------------------------------------
    /** @brief Integration of splines (described by their spline coefficients).
     *
     * See @ref MultipatchSplineEvaluatorOperator.
     *
     * @param[out] integrals An Kokkos::View on host containing the integrals of each spline on each patch. 
     * @param[in] patches_splines A MultipatchType of DField storing the 2D spline coefficients.
     */
    void integrate(
            DKokkosView_h<n_patches> const& integrals,
            MultipatchSplineCoeff const& patches_splines) const
    {
        (apply_integrate<Patches>(
                 integrals(ddc::type_seq_rank_v<Patches, PatchOrdering>),
                 patches_splines.template get<Patches>()),
         ...);
    }


public:
    // Apply functions to manage the values on each patch. ---------------------------------------

    /** @brief Compute the values or the derivatives of a given patch at the coordinates
     * defined on the given patch.
     * Needed public for functions on GPU. 
     * @tparam EvalType1 Evaluation type: either eval_type or eval_deriv_type.
     * @tparam EvalType2 Evaluation type: either eval_type or eval_deriv_type.
     * @tparam StoringPatch Patch type where the given coordinates are stored. 
     *      They are not especially physically located on this patch.
     * @param[out] patch_values Field of values of the function or derivative. 
     * @param[in] patch_coords ConstField of coordinates defined on the StoringPatch and 
     *          where we want to evaluate the function or derivative. 
     * @param[in] patches_splines MultipatchType of spline coefficients of the splines 
     *          on every patches. 
     */
    template <class EvalType1, class EvalType2, class StoringPatch>
    void apply_evaluator(
            ValuesOnPatch<StoringPatch> const& patch_values,
            CoordConstFieldOnPatch<StoringPatch> const& patch_coords,
            MultipatchSplineCoeff const patches_splines) const
    {
        assert(get_idx_range(patch_values) == get_idx_range(patch_coords));

        using Index =
                typename batched_evaluation_idx_range_type<StoringPatch>::discrete_element_type;
        batched_evaluation_idx_range_type<StoringPatch> idx_range = get_idx_range(patch_values);

        ddc::parallel_for_each(
                exec_space(),
                idx_range,
                KOKKOS_CLASS_LAMBDA(Index const& idx) {
                    CoordOnPatch<StoringPatch> const coord = patch_coords(idx);
                    int const patch_idx = get_patch_idx(coord);
                    if (patch_idx < 0
                        && !((std::is_same_v<EvalType1, eval_type>)&&(
                                std::is_same_v<EvalType2, eval_type>))) {
                        Kokkos::abort("The evaluation coordinate has to be on a patch."
                                      "No extrapolation rule for derivatives. \n");
                    }
                    patch_values(idx) = recursive_dispatch_patch_function<
                            EvalType1,
                            EvalType2>(coord, patches_splines, patch_idx);
                });
    }

    /** @brief Integrate the spline defined on the given patch.
     * @tparam Patch Patch type where the integration of the spline is computed. 
     * @param[out] integral Double, value of the integral of the spline on the given Patch.
     * @param[in] spline_coef ConstField of spline coefficients of the spline on the given Patch.
     */
    template <class Patch>
    void apply_integrate(double& integral, SplineCoeffOnPatch<Patch> const& spline_coef) const
    {
        static_assert(
                std::is_same_v<exec_space, Kokkos::DefaultHostExecutionSpace>,
                "Can only be called on CPU: .integrals() not defined yet on GPU.");
        static_assert(
                Kokkos::SpaceAccessibility<exec_space, memory_space>::accessible,
                "Execution space and memory space have to be compatible.");

        using bsplines_1 = bsplines_type1<Patch>;
        using bsplines_2 = bsplines_type2<Patch>;
        using IdxBS12 = typename spline_idx_range_type<Patch>::discrete_element_type;

        DFieldMem<IdxRange<bsplines_1>, memory_space> values1_alloc(
                get_idx_range<bsplines_1>(spline_coef));
        DField<IdxRange<bsplines_1>, memory_space> values1 = get_field(values1_alloc);
        DFieldMem<IdxRange<bsplines_2>, memory_space> values2_alloc(
                get_idx_range<bsplines_2>(spline_coef));
        DField<IdxRange<bsplines_2>, memory_space> values2 = get_field(values2_alloc);
        ddc::integrals(exec_space(), values1);
        ddc::integrals(exec_space(), values2);

        integral = ddc::parallel_transform_reduce(
                exec_space(),
                get_idx_range(spline_coef),
                0.0,
                ddc::reducer::sum<double>(),
                KOKKOS_LAMBDA(IdxBS12 const i12) {
                    return spline_coef(i12) * values1(ddc::select<bsplines_1>(i12))
                           * values2(ddc::select<bsplines_2>(i12));
                });
    }



private:
    // Recursive method to dispatch the coordinates on the right patch ---------------------------

    /// @brief Dispatch the given coordinate on the right patch to evaluate the right spline.
    template <class EvalType1, class EvalType2, class Dim1, class Dim2, int TestPatchIdx = 0>
    KOKKOS_INLINE_FUNCTION double recursive_dispatch_patch_function(
            Coord<Dim1, Dim2> coord,
            MultipatchSplineCoeff const& patches_splines,
            int const patch_idx) const
    {
        if constexpr (TestPatchIdx == ddc::type_seq_size_v<PatchOrdering>) {
            if (patch_idx >= 0) {
                Kokkos::abort("The recursion has reached the end without finding where "
                              "the coordinate is physically located.");
            }
            // Coord not on patch. Stop recursing.
            if constexpr (
                    std::is_same_v<EvalType1, eval_type> && std::is_same_v<EvalType2, eval_type>) {
                /* The operator currently works only for the case where the continuous 
                   dimensions of all the patches are the same. So the equivalent coordinates
                   are the same on any patches. 
                   Add a method to determine the patch from continuous dimensions, once 
                   we would like to deal with other cases. 
                */
                using AnyPatch = ddc::type_seq_element_t<0, PatchOrdering>;
                replace_periodic_coord_inside<AnyPatch>(coord);
                return m_extrapolation_rule(coord, patches_splines, patch_idx);
            } else {
                Kokkos::abort("The spline derivatives cannot be evaluated at coordinates "
                              "outside of the domain.");
            }
        } else {
            if (patch_idx == TestPatchIdx) {
                using TestPatch = ddc::type_seq_element_t<TestPatchIdx, PatchOrdering>;
                CoordOnPatch<TestPatch> test_coord = get_equivalent_coord<
                        typename TestPatch::Dim1,
                        typename TestPatch::Dim2,
                        Dim1,
                        Dim2>(coord);
                replace_periodic_coord_inside<TestPatch>(test_coord);

                SplineCoeffOnPatch<TestPatch> const test_spline
                        = patches_splines.template get<TestPatch>();
                return eval_no_bc<EvalType1, EvalType2, TestPatch>(test_coord, test_spline);
            } else {
                return recursive_dispatch_patch_function<
                        EvalType1,
                        EvalType2,
                        Dim1,
                        Dim2,
                        TestPatchIdx + 1>(coord, patches_splines, patch_idx);
            }
        }
    }


    // Useful functions --------------------------------------------------------------------------

    /// @brief Call the patch locator to get the index of the patch where the given coordinate
    /// is physically located.
    template <class Dim1, class Dim2>
    KOKKOS_INLINE_FUNCTION int get_patch_idx(Coord<Dim1, Dim2> const coord) const
    {
        using Mapping = typename PatchLocator::get_mapping_on_logical_dim_t<Dim1, Dim2>;
        Mapping const mapping(m_patch_locator.template get_mapping_on_logical_dim<Dim1, Dim2>());
        return m_patch_locator(mapping(coord));
    }

    /// @brief Call mappings to get the equivalent coordinate defined on a current patch
    /// on the target patch. Pass by the physical domain.
    template <class TargetDim1, class TargetDim2, class CurrentDim1, class CurrentDim2>
    KOKKOS_INLINE_FUNCTION Coord<TargetDim1, TargetDim2> get_equivalent_coord(
            Coord<CurrentDim1, CurrentDim2> const& current_coord) const
    {
        if constexpr (std::is_same_v<
                              Coord<TargetDim1, TargetDim2>,
                              Coord<CurrentDim1, CurrentDim2>>) {
            return current_coord;
        } else {
            using CurrentMapping =
                    typename PatchLocator::get_mapping_on_logical_dim_t<CurrentDim1, CurrentDim2>;
            using TargetMapping =
                    typename PatchLocator::get_mapping_on_logical_dim_t<TargetDim1, TargetDim2>;

            static_assert((std::is_same_v<
                           Coord<typename CurrentMapping::curvilinear_tag_r,
                                 typename CurrentMapping::curvilinear_tag_theta>,
                           Coord<CurrentDim1, CurrentDim2>>));
            static_assert((std::is_same_v<
                           Coord<typename TargetMapping::curvilinear_tag_r,
                                 typename TargetMapping::curvilinear_tag_theta>,
                           Coord<TargetDim1, TargetDim2>>));

            CurrentMapping const current_mapping(
                    m_patch_locator
                            .template get_mapping_on_logical_dim<CurrentDim1, CurrentDim2>());
            TargetMapping const target_mapping(
                    m_patch_locator.template get_mapping_on_logical_dim<TargetDim1, TargetDim2>());

            return target_mapping(current_mapping(current_coord));
        }
    }

    /// @brief Replace a coordinate inside the domain if it is periodic.
    template <class Patch>
    KOKKOS_INLINE_FUNCTION void replace_periodic_coord_inside(CoordOnPatch<Patch>& coord) const
    {
        using bsplines_1 = bsplines_type1<Patch>;
        using bsplines_2 = bsplines_type2<Patch>;

        using Dim1 = continuous_dimension_type1<Patch>;
        using Dim2 = continuous_dimension_type2<Patch>;

        Coord<Dim1> coord_1(coord);
        Coord<Dim2> coord_2(coord);

        ddcHelper::restrict_to_bspline_domain<bsplines_1>(coord_1);
        ddcHelper::restrict_to_bspline_domain<bsplines_2>(coord_2);

        coord = CoordOnPatch<Patch>(coord_1, coord_2);
    }


    // Evaluation functions ----------------------------------------------------------------------

    /// @brief Evaluate the given spline at the given coordinate without carring of the boundary
    /// conditions.
    template <class EvalType1, class EvalType2, class Patch, class Layout>
    KOKKOS_INLINE_FUNCTION double eval_no_bc(
            CoordOnPatch<Patch> const& coord_eval,
            DConstField<spline_idx_range_type<Patch>, memory_space, Layout> const& spline_coef)
            const
    {
        static_assert(
                std::is_same_v<EvalType1, eval_type> || std::is_same_v<EvalType1, eval_deriv_type>);
        static_assert(
                std::is_same_v<EvalType2, eval_type> || std::is_same_v<EvalType2, eval_deriv_type>);

        using bsplines_1 = bsplines_type1<Patch>;
        using bsplines_2 = bsplines_type2<Patch>;

        Idx<bsplines_1> jmin1;
        Idx<bsplines_2> jmin2;

        std::array<double, bsplines_1::degree() + 1> vals1_ptr;
        DSpan1D const vals1(vals1_ptr.data(), bsplines_1::degree() + 1);
        std::array<double, bsplines_2::degree() + 1> vals2_ptr;
        DSpan1D const vals2(vals2_ptr.data(), bsplines_2::degree() + 1);

        Coord<continuous_dimension_type1<Patch>> coord_eval_interest1
                = ddc::select<continuous_dimension_type1<Patch>>(coord_eval);
        Coord<continuous_dimension_type2<Patch>> coord_eval_interest2
                = ddc::select<continuous_dimension_type2<Patch>>(coord_eval);

        if constexpr (std::is_same_v<EvalType1, eval_type>) {
            jmin1 = ddc::discrete_space<bsplines_1>().eval_basis(vals1, coord_eval_interest1);
        } else if constexpr (std::is_same_v<EvalType1, eval_deriv_type>) {
            jmin1 = ddc::discrete_space<bsplines_1>().eval_deriv(vals1, coord_eval_interest1);
        }

        if constexpr (std::is_same_v<EvalType2, eval_type>) {
            jmin2 = ddc::discrete_space<bsplines_2>().eval_basis(vals2, coord_eval_interest2);
        } else if constexpr (std::is_same_v<EvalType2, eval_deriv_type>) {
            jmin2 = ddc::discrete_space<bsplines_2>().eval_deriv(vals2, coord_eval_interest2);
        }

        double y = 0.0;
        for (std::size_t i = 0; i < bsplines_1::degree() + 1; ++i) {
            for (std::size_t j = 0; j < bsplines_2::degree() + 1; ++j) {
                y += spline_coef(jmin1 + i, jmin2 + j) * vals1[i] * vals2[j];
            }
        }
        return y;
    }
};
