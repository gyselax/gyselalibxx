

# File multipatch\_spline\_builder\_2d.hpp

[**File List**](files.md) **>** [**multipatch**](dir_7740c6927b2da0a836b00bedb040a06d.md) **>** [**spline**](dir_729d943c83b6b5573a69e28a4db4673a.md) **>** [**multipatch\_spline\_builder\_2d.hpp**](multipatch__spline__builder__2d_8hpp.md)

[Go to the documentation of this file](multipatch__spline__builder__2d_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>
#include <tuple>
#include <utility>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"
#include "multipatch_type.hpp"
#include "spline_builder_deriv_field_2d.hpp"
#include "types.hpp"

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
        ddc::SplineBuilderClosure SBCLower1,
        ddc::SplineBuilderClosure SBCUpper1,
        ddc::SplineBuilderClosure SBCLower2,
        ddc::SplineBuilderClosure SBCUpper2,
        ddc::SplineBuilderClosure BcTransition,
        class Connectivity,
        ddc::SplineSolver Solver,
        template <typename P>
        typename ValuesOnPatch,
        class... Patches>
class MultipatchSplineBuilder2D
{
    static_assert(
            (((ddc::is_uniform_bsplines_v<BSpline1OnPatch<Patches>>)
              || (ddc::is_non_uniform_bsplines_v<BSpline1OnPatch<Patches>>))
             && ...),
            "The BSpline1OnPatch argument does not create 1D BSpline objects.");
    static_assert(
            (((ddc::is_uniform_bsplines_v<BSpline2OnPatch<Patches>>)
              || (ddc::is_non_uniform_bsplines_v<BSpline2OnPatch<Patches>>))
             && ...),
            "The BSpline2OnPatch argument does not create 1D BSpline objects.");
    static_assert(
            (((ddc::is_uniform_point_sampling_v<Grid1OnPatch<Patches>>)
              || (ddc::is_non_uniform_point_sampling_v<Grid1OnPatch<Patches>>))
             && ...),
            "The Grid1OnPatch argument does not create 1D Grid objects.");
    static_assert(
            (((ddc::is_uniform_point_sampling_v<Grid2OnPatch<Patches>>)
              || (ddc::is_non_uniform_point_sampling_v<Grid2OnPatch<Patches>>))
             && ...),
            "The Grid2OnPatch argument does not create 1D Grid objects.");
    static_assert(
            ((std::is_same_v<
                    typename BSpline1OnPatch<Patches>::continuous_dimension_type,
                    typename Grid1OnPatch<Patches>::continuous_dimension_type>)&&...),
            "The BSpline1OnPatch argument does not define B-splines on the dimension where the "
            "grids defined by Grid1OnPatch are defined.");
    static_assert(
            ((std::is_same_v<
                    typename BSpline2OnPatch<Patches>::continuous_dimension_type,
                    typename Grid2OnPatch<Patches>::continuous_dimension_type>)&&...),
            "The BSpline2OnPatch argument does not define B-splines on the dimension where the "
            "grids defined by Grid2OnPatch are defined.");
    static_assert(
            !(is_deriv_field_v<ValuesOnPatch<Patches>> && ...),
            "ValuesOnPatch represents the field, not the field with derivatives.");

    template <class Patch, class FieldType>
    struct Build_BuilderType
    {
        static_assert(
                !std::is_same_v<Patch, Patch>,
                "The values should be saved in a constant field of doubles on the specified memory "
                "space.");
    };

    template <class Patch, class... Grid1D>
    struct Build_BuilderType<Patch, DConstField<IdxRange<Grid1D...>, MemorySpace>>
    {
        using lower_matching_edge1 = equivalent_edge_t<
                Edge<Patch, Grid1OnPatch<Patch>, FRONT>,
                typename Connectivity::interface_collection>;
        using upper_matching_edge1 = equivalent_edge_t<
                Edge<Patch, Grid1OnPatch<Patch>, BACK>,
                typename Connectivity::interface_collection>;
        using lower_matching_edge2 = equivalent_edge_t<
                Edge<Patch, Grid2OnPatch<Patch>, FRONT>,
                typename Connectivity::interface_collection>;
        using upper_matching_edge2 = equivalent_edge_t<
                Edge<Patch, Grid2OnPatch<Patch>, BACK>,
                typename Connectivity::interface_collection>;
        using type = ddc::SplineBuilder2D<
                ExecSpace,
                MemorySpace,
                BSpline1OnPatch<Patch>,
                BSpline2OnPatch<Patch>,
                Grid1OnPatch<Patch>,
                Grid2OnPatch<Patch>,
                std::is_same_v<lower_matching_edge1, OutsideEdge> ? SBCLower1 : BcTransition,
                std::is_same_v<upper_matching_edge1, OutsideEdge> ? SBCUpper1 : BcTransition,
                std::is_same_v<lower_matching_edge2, OutsideEdge> ? SBCLower2 : BcTransition,
                std::is_same_v<upper_matching_edge2, OutsideEdge> ? SBCUpper2 : BcTransition,
                Solver>;
    };

    template <class Patch>
    using BuilderOnPatch = typename Build_BuilderType<Patch, ValuesOnPatch<Patch>>::type;

    template <class Patch>
    using SplineOnPatch = DField<
            typename BuilderOnPatch<Patch>::template batched_spline_domain_type<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    template <class Patch>
    using Derivs1OnPatch = DConstField<
            typename BuilderOnPatch<Patch>::template batched_derivs_domain_type1<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    template <class Patch>
    using Derivs2OnPatch = DConstField<
            typename BuilderOnPatch<Patch>::template batched_derivs_domain_type2<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    template <class Patch>
    using Derivs12OnPatch = DConstField<
            typename BuilderOnPatch<Patch>::template batched_derivs_domain_type<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    template <class Patch>
    using DerivConstFieldOnPatch = DerivConstField<
            double,
            IdxRange<
                    ddc::Deriv<typename Patch::Dim1>,
                    typename Patch::Grid1,
                    ddc::Deriv<typename Patch::Dim2>,
                    typename Patch::Grid2>,
            MemorySpace>;

    template <class Patch, class BatchedIdxRange>
    using batched_spline_domain_type_on_patch
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<BatchedIdxRange>,
                    ddc::detail::TypeSeq<typename Patch::Grid1, typename Patch::Grid2>,
                    ddc::detail::TypeSeq<typename Patch::BSplines1, typename Patch::BSplines2>>>;

    template <class Patch, class BatchedIdxRange>
    using batched_derivs_domain_type_on_patch
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<BatchedIdxRange>,
                    ddc::detail::TypeSeq<typename Patch::Grid1, typename Patch::Grid2>,
                    ddc::detail::TypeSeq<
                            ddc::Deriv<typename Patch::Dim1>,
                            ddc::Deriv<typename Patch::Dim2>>>>;


    using MultipatchSplineCoeffs = MultipatchField<SplineOnPatch, Patches...>;

    using MultipatchValues = MultipatchField<ValuesOnPatch, Patches...>;

    using MultipatchDerivs1 = MultipatchField<Derivs1OnPatch, Patches...>;

    using MultipatchDerivs2 = MultipatchField<Derivs2OnPatch, Patches...>;

    using MultipatchDerivs12 = MultipatchField<Derivs12OnPatch, Patches...>;

    using MultipatchDerivField = MultipatchField<DerivConstFieldOnPatch, Patches...>;

    using MultipatchIdxRange = MultipatchType<IdxRangeOnPatch, Patches...>;

    using MultipatchIdxRangeBS = MultipatchType<IdxRangeBSOnPatch, Patches...>;

    template <class MultipatchIdxRangeBatched>
    using MultipatchIdxRangeBatchedBS = detail::MultipatchType<
            typename MultipatchIdxRangeBatched::PatchOrdering,
            ddc::detail::TypeSeq<batched_spline_domain_type_on_patch<
                    Patches,
                    typename MultipatchIdxRangeBatched::template TypeOnPatch<Patches>>...>>;

    template <class MultipatchIdxRangeBatched>
    using MultipatchIdxRangeBatchedDeriv = detail::MultipatchType<
            typename MultipatchIdxRangeBatched::PatchOrdering,
            ddc::detail::TypeSeq<batched_derivs_domain_type_on_patch<
                    Patches,
                    typename MultipatchIdxRangeBatched::template TypeOnPatch<Patches>>...>>;

    using BuilderTuple = std::tuple<BuilderOnPatch<Patches> const&...>;


    BuilderTuple const m_builders;

public:
    using exec_space = ExecSpace;

    using memory_space = MemorySpace;

    using data_type = double;

    using interpolation_idx_range_type = MultipatchIdxRange;

    using coeff_idx_range_type = MultipatchIdxRangeBS;

    static constexpr std::size_t rank()
    {
        return 2;
    }

    template <class MultipatchIdxRangeBatched>
    using batched_basis_idx_range_type = MultipatchIdxRangeBatchedBS<MultipatchIdxRangeBatched>;

    template <class MultipatchIdxRangeBatched>
    using batched_derivs_idx_range_type = MultipatchIdxRangeBatchedDeriv<MultipatchIdxRangeBatched>;

private:
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
    explicit MultipatchSplineBuilder2D(BuilderOnPatch<Patches> const&... builders)
        : m_builders(std::tie(builders...)) {};


    ~MultipatchSplineBuilder2D() = default;

    void operator()(
            MultipatchSplineCoeffs splines,
            MultipatchValues const& values,
            std::optional<MultipatchDerivs1> derivs_min1 = std::nullopt,
            std::optional<MultipatchDerivs1> derivs_max1 = std::nullopt,
            std::optional<MultipatchDerivs2> derivs_min2 = std::nullopt,
            std::optional<MultipatchDerivs2> derivs_max2 = std::nullopt,
            std::optional<MultipatchDerivs12> cross_derivs_min1_min2 = std::nullopt,
            std::optional<MultipatchDerivs12> cross_derivs_max1_min2 = std::nullopt,
            std::optional<MultipatchDerivs12> cross_derivs_min1_max2 = std::nullopt,
            std::optional<MultipatchDerivs12> cross_derivs_max1_max2 = std::nullopt) const
    {
        ((std::get<BuilderOnPatch<Patches> const&>(m_builders)(
                 splines.template get<Patches>(),
                 get_const_field(values.template get<Patches>()),
                 get_deriv_value<Patches>(derivs_min1),
                 get_deriv_value<Patches>(derivs_max1),
                 get_deriv_value<Patches>(derivs_min2),
                 get_deriv_value<Patches>(derivs_max2),
                 get_deriv_value<Patches>(cross_derivs_min1_min2),
                 get_deriv_value<Patches>(cross_derivs_max1_min2),
                 get_deriv_value<Patches>(cross_derivs_min1_max2),
                 get_deriv_value<Patches>(cross_derivs_max1_max2))),
         ...);
    };

    void operator()(MultipatchSplineCoeffs splines, MultipatchDerivField functions_and_derivs) const
    {
        (apply_builder<Patches>(
                 std::get<BuilderOnPatch<Patches> const&>(m_builders),
                 splines.template get<Patches>(),
                 functions_and_derivs.template get<Patches>()),
         ...);
    };

    MultipatchIdxRangeBS const spline_idx_ranges() const
    {
        return MultipatchIdxRangeBS(
                get_spline_idx_range(std::get<BuilderOnPatch<Patches> const&>(m_builders))...);
    };

    template <class MultipatchIdxRangeBatched>
    MultipatchIdxRangeBatchedBS<MultipatchIdxRangeBatched> const batched_basis_idx_ranges(
            MultipatchIdxRangeBatched idx_ranges) const
    {
        return MultipatchIdxRangeBatchedBS<MultipatchIdxRangeBatched>(batched_basis_idx_range(
                std::get<BuilderOnPatch<Patches> const&>(m_builders),
                idx_ranges.template get<Patches>())...);
    };

private:
    template <class PatchP>
    static void apply_builder(
            BuilderOnPatch<PatchP> const& builder,
            SplineOnPatch<PatchP> spline,
            DerivConstFieldOnPatch<PatchP> function_and_deriv)
    {
        SplineBuilderDerivField2D<
                ExecSpace,
                BSpline1OnPatch<PatchP>,
                BSpline2OnPatch<PatchP>,
                Grid1OnPatch<PatchP>,
                Grid2OnPatch<PatchP>,
                BuilderOnPatch<PatchP>::builder_type1::s_sbc_xmin,
                BuilderOnPatch<PatchP>::builder_type1::s_sbc_xmax,
                BuilderOnPatch<PatchP>::builder_type2::s_sbc_xmin,
                BuilderOnPatch<PatchP>::builder_type2::s_sbc_xmax>
                builder_applier(builder);
        builder_applier(spline, function_and_deriv);
    };
};
```


