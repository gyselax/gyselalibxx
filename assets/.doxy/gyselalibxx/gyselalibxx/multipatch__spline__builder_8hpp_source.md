

# File multipatch\_spline\_builder.hpp

[**File List**](files.md) **>** [**multipatch**](dir_7740c6927b2da0a836b00bedb040a06d.md) **>** [**spline**](dir_729d943c83b6b5573a69e28a4db4673a.md) **>** [**multipatch\_spline\_builder.hpp**](multipatch__spline__builder_8hpp.md)

[Go to the documentation of this file](multipatch__spline__builder_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>
#include <tuple>
#include <utility>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"

template <
        class ExecSpace,
        class MemorySpace,
        template <typename P>
        typename BSplineOnPatch,
        template <typename P>
        typename GridOnPatch,
        ddc::BoundCond BcLower,
        ddc::BoundCond BcUpper,
        ddc::BoundCond BcTransition,
        class Connectivity,
        ddc::SplineSolver Solver,
        template <typename P>
        typename ValuesOnPatch,
        class... Patches>
class MultipatchSplineBuilder
{
    static_assert(
            (((ddc::is_uniform_bsplines_v<BSplineOnPatch<Patches>>)
              || (ddc::is_non_uniform_bsplines_v<BSplineOnPatch<Patches>>))
             && ...),
            "The BSplineOnPatch argument does not create 1D BSpline objects.");
    static_assert(
            (((ddc::is_uniform_point_sampling_v<GridOnPatch<Patches>>)
              || (ddc::is_non_uniform_point_sampling_v<GridOnPatch<Patches>>))
             && ...),
            "The GridOnPatch argument does not create 1D Grid objects.");

    using PatchOrdering = ddc::detail::TypeSeq<Patches...>;
    static constexpr std::size_t n_patches = ddc::type_seq_size_v<PatchOrdering>;

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
        using lower_matching_edge = equivalent_edge_t<
                Edge<Patch, GridOnPatch<Patch>, FRONT>,
                typename Connectivity::interface_collection>;
        using upper_matching_edge = equivalent_edge_t<
                Edge<Patch, GridOnPatch<Patch>, BACK>,
                typename Connectivity::interface_collection>;
        static constexpr std::size_t patch_id = ddc::type_seq_rank_v<Patch, PatchOrdering>;
        using type = ddc::SplineBuilder<
                ExecSpace,
                MemorySpace,
                BSplineOnPatch<Patch>,
                GridOnPatch<Patch>,
                std::is_same_v<lower_matching_edge, OutsideEdge> ? BcLower : BcTransition,
                std::is_same_v<upper_matching_edge, OutsideEdge> ? BcUpper : BcTransition,
                Solver,
                Grid1D...>;
    };

    template <class Patch>
    using BuilderOnPatch = typename Build_BuilderType<Patch, ValuesOnPatch<Patch>>::type;

    template <class Patch>
    using SplineOnPatch
            = DField<typename BuilderOnPatch<Patch>::batched_spline_domain_type, MemorySpace>;

    template <class Patch>
    using DerivsOnPatch
            = DConstField<typename BuilderOnPatch<Patch>::batched_derivs_domain_type, MemorySpace>;

    using MultipatchSplineCoeffs = MultipatchField<SplineOnPatch, Patches...>;

    // For PERIODIC or GREVILLE boundary conditions
    using MultipatchValues = MultipatchField<ValuesOnPatch, Patches...>;

    using MultipatchDerivs = MultipatchField<DerivsOnPatch, Patches...>;

    using BuilderTuple = std::tuple<BuilderOnPatch<Patches> const&...>;


    BuilderTuple const m_builders;

private:
    template <class Patch>
    std::optional<DerivsOnPatch<Patch>> get_deriv_value(
            std::optional<MultipatchDerivs> derivs) const
    {
        if (derivs.has_value()) {
            return derivs->template get<Patch>();
        } else {
            return std::nullopt;
        }
    }

public:
    explicit MultipatchSplineBuilder(BuilderOnPatch<Patches> const&... builders)
        : m_builders(std::tie(builders...)) {};


    ~MultipatchSplineBuilder() = default;

    void operator()(
            MultipatchSplineCoeffs splines,
            MultipatchValues const& values,
            std::optional<MultipatchDerivs> derivs_xmin = std::nullopt,
            std::optional<MultipatchDerivs> derivs_xmax = std::nullopt) const
    {
        ((std::get<BuilderOnPatch<Patches> const&>(m_builders)(
                 splines.template get<Patches>(),
                 get_const_field(values.template get<Patches>()),
                 get_deriv_value<Patches>(derivs_xmin),
                 get_deriv_value<Patches>(derivs_xmax))),
         ...);
    };
};
```


