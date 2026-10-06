

# File multipatch\_spline\_evaluator\_2d.hpp

[**File List**](files.md) **>** [**multipatch**](dir_7740c6927b2da0a836b00bedb040a06d.md) **>** [**spline**](dir_729d943c83b6b5573a69e28a4db4673a.md) **>** [**multipatch\_spline\_evaluator\_2d.hpp**](multipatch__spline__evaluator__2d_8hpp.md)

[Go to the documentation of this file](multipatch__spline__evaluator__2d_8hpp.md)


```C++
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
#include "view.hpp"


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
    template <class Patch>
    using continuous_dimension_type1 = typename BSpline1OnPatch<Patch>::continuous_dimension_type;

    template <class Patch>
    using continuous_dimension_type2 = typename BSpline2OnPatch<Patch>::continuous_dimension_type;


    using exec_space = ExecSpace;

    using memory_space = MemorySpace;


    template <class Patch>
    using evaluation_discrete_dimension_type1 = Grid1OnPatch<Patch>;

    template <class Patch>
    using evaluation_discrete_dimension_type2 = Grid2OnPatch<Patch>;

    template <class Patch>
    using bsplines_type1 = BSpline1OnPatch<Patch>;

    template <class Patch>
    using bsplines_type2 = BSpline2OnPatch<Patch>;


    template <class Patch>
    using evaluation_idx_range_type1 = IdxRange<evaluation_discrete_dimension_type1<Patch>>;

    template <class Patch>
    using evaluation_idx_range_type2 = IdxRange<evaluation_discrete_dimension_type2<Patch>>;

    template <class Patch>
    using evaluation_idx_range_type = IdxRange<
            evaluation_discrete_dimension_type1<Patch>,
            evaluation_discrete_dimension_type2<Patch>>;


    template <class Patch>
    using batched_evaluation_idx_range_type = typename ValuesOnPatch<Patch>::discrete_domain_type;


    template <class Patch>
    using spline_idx_range_type1 = IdxRange<bsplines_type1<Patch>>;

    template <class Patch>
    using spline_idx_range_type2 = IdxRange<bsplines_type2<Patch>>;

    template <class Patch>
    using spline_idx_range_type = IdxRange<bsplines_type1<Patch>, bsplines_type2<Patch>>;

public:
    template <class Patch>
    using SplineCoeffOnPatch = DConstField<spline_idx_range_type<Patch>, MemorySpace>;

private:
    // MultipatchTypes
    using MultipatchValues = MultipatchField<ValuesOnPatch, Patches...>;
    using MultipatchCoordField = MultipatchField<CoordConstFieldOnPatch, Patches...>;
    using MultipatchSplineCoeff = MultipatchField<SplineCoeffOnPatch, Patches...>;

    // Patches
    using PatchOrdering = ddc::detail::TypeSeq<Patches...>;


    // As we work with the same continuous dimensions --------------------------------------------
    // Continuous first dimension of the patches.
    using CDim1 = continuous_dimension_type1<ddc::type_seq_element_t<0, PatchOrdering>>;
    // Continuous second dimension of the patches.
    using CDim2 = continuous_dimension_type2<ddc::type_seq_element_t<0, PatchOrdering>>;

    // Dimension for the derivatives along the first dimension.
    using DerivDim1 = ddc::Deriv<CDim1>;
    // Dimension for the derivatives along the second dimension.
    using DerivDim2 = ddc::Deriv<CDim2>;

private:
    // Fields
    template <class Patch>
    using CoordConstFieldOnPatch
            = ConstField<Coord<CDim1, CDim2>, evaluation_idx_range_type<Patch>, MemorySpace>;

public:
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
            "continuous dimension. E.g. OnionPatchLocator.");


    // Members -----------------------------------------------------------------------------------
    PatchLocator const m_patch_locator;
    ExtrapolationRule const m_extrapolation_rule;


public:
    ~MultipatchSplineEvaluator2D() = default;

    MultipatchSplineEvaluator2D(
            PatchLocator const& patch_locator,
            ExtrapolationRule const& extrapolation_rule)
        : m_patch_locator(patch_locator)
        , m_extrapolation_rule(extrapolation_rule)
    {
    }


    // Assignment operators ----------------------------------------------------------------------
    template <class Coord>
    KOKKOS_FUNCTION double operator()(
            Coord const coord_eval,
            MultipatchSplineCoeff const& patches_splines) const
    {
        int const patch_idx = get_patch_idx(coord_eval);
        return recursive_dispatch_patch_function(Idx<>(), coord_eval, patches_splines, patch_idx);
    }

    void operator()(
            MultipatchValues const& patches_values,
            MultipatchCoordField const& patches_coords,
            MultipatchSplineCoeff const& patches_splines) const
    {
        (apply_evaluator<Patches>(
                 Idx<>(),
                 patches_values.template get<Patches>(),
                 patches_coords.template get<Patches>(),
                 patches_splines),
         ...);
    }


    // Derivatives operators ---------------------------------------------------------------------
    template <class... DerivDims>
    KOKKOS_FUNCTION double deriv(
            Idx<DerivDims...> const& deriv_order,
            Coord<CDim1, CDim2> const& coord_eval,
            MultipatchSplineCoeff const& patches_splines) const
    {
        static_assert(
                (ddc::in_tags_v<DerivDims, ddc::detail::TypeSeq<DerivDim1, DerivDim2>> && ...),
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
                patches_splines,
                patch_idx);
    }

    template <class... DerivDims>
    void deriv(
            Idx<DerivDims...> const& deriv_order,
            MultipatchValues const& patches_deriv,
            MultipatchCoordField const& patches_coords,
            MultipatchSplineCoeff const& patches_splines) const
    {
        static_assert(
                (ddc::in_tags_v<DerivDims, ddc::detail::TypeSeq<DerivDim1, DerivDim2>> && ...),
                "Please provide a deriv_order indexed on the derivative dimensions of the patches "
                "(same dimensions for all the patches).");
        (apply_evaluator<Patches>(
                 deriv_order,
                 patches_deriv.template get<Patches>(),
                 patches_coords.template get<Patches>(),
                 patches_splines),
         ...);
    }

    // Integrate operator ------------------------------------------------------------------------
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

    template <class StoringPatch, class IdxDeriv>
    void apply_evaluator(
            IdxDeriv const& deriv_order,
            ValuesOnPatch<StoringPatch> const& patch_values,
            CoordConstFieldOnPatch<StoringPatch> const& patch_coords,
            MultipatchSplineCoeff const patches_splines) const
    {
        assert(get_idx_range(patch_values) == get_idx_range(patch_coords));

        using Index =
                typename batched_evaluation_idx_range_type<StoringPatch>::discrete_element_type;
        batched_evaluation_idx_range_type<StoringPatch> idx_range = get_idx_range(patch_values);

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space(),
                idx_range,
                KOKKOS_CLASS_LAMBDA(Index const& idx) {
                    Coord<CDim1, CDim2> const coord = patch_coords(idx);
                    int const patch_idx = get_patch_idx(coord);
                    if (patch_idx < 0 && !(std::is_same_v<IdxDeriv, Idx<>>)) {
                        Kokkos::abort("The evaluation coordinate has to be on a patch."
                                      "No extrapolation rule for derivatives. \n");
                    }
                    patch_values(idx) = recursive_dispatch_patch_function(
                            deriv_order,
                            coord,
                            patches_splines,
                            patch_idx);
                });
    }

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
                "values1 (MultipatchSplineEvaluator2D::apply_integrate)",
                get_idx_range<bsplines_1>(spline_coef));
        DField<IdxRange<bsplines_1>, memory_space> values1 = get_field(values1_alloc);
        DFieldMem<IdxRange<bsplines_2>, memory_space> values2_alloc(
                "values2 (MultipatchSplineEvaluator2D::apply_integrate)",
                get_idx_range<bsplines_2>(spline_coef));
        DField<IdxRange<bsplines_2>, memory_space> values2 = get_field(values2_alloc);
        ddc::integrals(exec_space(), values1);
        ddc::integrals(exec_space(), values2);

        const std::source_location location = std::source_location::current();
        integral = ddc::parallel_transform_reduce(
                location.function_name(),
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

    template <class IdxDeriv, int TestPatchIdx = 0>
    KOKKOS_INLINE_FUNCTION double recursive_dispatch_patch_function(
            IdxDeriv const& deriv_order,
            Coord<CDim1, CDim2> coord,
            MultipatchSplineCoeff const& patches_splines,
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
                return m_extrapolation_rule(coord, patches_splines, patch_idx);
            } else {
                Kokkos::abort( // cppcheck-suppress missingReturn
                        "The spline derivatives cannot be evaluated at coordinates "
                        "outside of the domain.");
            }
        } else {
            if (patch_idx == TestPatchIdx) {
                using TestPatch = ddc::type_seq_element_t<TestPatchIdx, PatchOrdering>;
                replace_periodic_coord_inside<TestPatch>(coord);

                SplineCoeffOnPatch<TestPatch> const test_spline
                        = patches_splines.template get<TestPatch>();
                return eval_no_bc<TestPatch>(deriv_order, coord, test_spline);
            } else {
                return recursive_dispatch_patch_function<
                        IdxDeriv,
                        TestPatchIdx + 1>(deriv_order, coord, patches_splines, patch_idx);
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
        using bsplines_1 = bsplines_type1<Patch>;
        using bsplines_2 = bsplines_type2<Patch>;

        Coord<CDim1> coord_1(coord);
        Coord<CDim2> coord_2(coord);

        ddcHelper::restrict_to_bspline_domain<bsplines_1>(coord_1);
        ddcHelper::restrict_to_bspline_domain<bsplines_2>(coord_2);

        coord = Coord<CDim1, CDim2>(coord_1, coord_2);
    }


    // Evaluation functions ----------------------------------------------------------------------

    template <class Patch, class Layout, class... DerivDims>
    KOKKOS_INLINE_FUNCTION double eval_no_bc(
            Idx<DerivDims...> const& deriv_order,
            Coord<CDim1, CDim2> const& coord_eval,
            DConstField<spline_idx_range_type<Patch>, memory_space, Layout> const& spline_coef)
            const
    {
        using deriv_dims_seq = ddc::detail::TypeSeq<DerivDims...>;

        using bsplines_1 = bsplines_type1<Patch>;
        using bsplines_2 = bsplines_type2<Patch>;

        Idx<bsplines_1> jmin1;
        Idx<bsplines_2> jmin2;

        std::array<double, bsplines_1::degree() + 1> vals1_ptr;
        DSpan1D const vals1(vals1_ptr.data(), bsplines_1::degree() + 1);
        std::array<double, bsplines_2::degree() + 1> vals2_ptr;
        DSpan1D const vals2(vals2_ptr.data(), bsplines_2::degree() + 1);

        Coord<CDim1> coord_eval_interest1 = ddc::select<CDim1>(coord_eval);
        Coord<CDim2> coord_eval_interest2 = ddc::select<CDim2>(coord_eval);

        if constexpr (!ddc::in_tags_v<DerivDim1, deriv_dims_seq>) {
            jmin1 = ddc::discrete_space<bsplines_1>().eval_basis(vals1, coord_eval_interest1);
        } else {
            const std::size_t order1 = (Idx<DerivDim1>(deriv_order) - Idx<DerivDim1>(0)).value();
            KOKKOS_ASSERT(order1 <= bsplines_1::degree())

            std::array<double, (bsplines_1::degree() + 1) * (bsplines_1::degree() + 1)> derivs1_ptr;
            Kokkos::mdspan<
                    double,
                    Kokkos::extents<
                            std::size_t,
                            bsplines_1::degree() + 1,
                            Kokkos::dynamic_extent>> const derivs1(derivs1_ptr.data(), order1 + 1);

            jmin1 = ddc::discrete_space<bsplines_1>()
                            .eval_basis_and_n_derivs(derivs1, coord_eval_interest1, order1);

            for (std::size_t i = 0; i < bsplines_1::degree() + 1; ++i) {
                vals1[i] = derivs1(i, order1);
            }
        }

        if constexpr (!ddc::in_tags_v<DerivDim2, deriv_dims_seq>) {
            jmin2 = ddc::discrete_space<bsplines_2>().eval_basis(vals2, coord_eval_interest2);
        } else {
            const std::size_t order2 = (Idx<DerivDim2>(deriv_order) - Idx<DerivDim2>(0)).value();
            KOKKOS_ASSERT(order2 <= bsplines_2::degree())

            std::array<double, (bsplines_2::degree() + 1) * (bsplines_2::degree() + 1)> derivs2_ptr;
            Kokkos::mdspan<
                    double,
                    Kokkos::extents<
                            std::size_t,
                            bsplines_2::degree() + 1,
                            Kokkos::dynamic_extent>> const derivs2(derivs2_ptr.data(), order2 + 1);

            jmin2 = ddc::discrete_space<bsplines_2>()
                            .eval_basis_and_n_derivs(derivs2, coord_eval_interest2, order2);

            for (std::size_t i = 0; i < bsplines_2::degree() + 1; ++i) {
                vals2[i] = derivs2(i, order2);
            }
        }


        double y = 0.0;
        for (std::size_t i = 0; i < bsplines_1::degree() + 1; ++i) {
            for (std::size_t j = 0; j < bsplines_2::degree() + 1; ++j) {
                Idx<bsplines_1, bsplines_2> idx_spline(jmin1 + i, jmin2 + j);
                y += spline_coef(idx_spline) * vals1[i] * vals2[j];
            }
        }
        return y;
    }
};
```


