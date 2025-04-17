

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
#include "multipatch_type.hpp"

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
        ddc::BoundCond BcLower1,
        ddc::BoundCond BcUpper1,
        ddc::BoundCond BcLower2,
        ddc::BoundCond BcUpper2,
        ddc::BoundCond BcTransition,
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
                std::is_same_v<lower_matching_edge1, OutsideEdge> ? BcLower1 : BcTransition,
                std::is_same_v<upper_matching_edge1, OutsideEdge> ? BcUpper1 : BcTransition,
                std::is_same_v<lower_matching_edge2, OutsideEdge> ? BcLower2 : BcTransition,
                std::is_same_v<upper_matching_edge2, OutsideEdge> ? BcUpper2 : BcTransition,
                Solver>;
    };

    template <class Patch>
    using BuilderOnPatch = typename Build_BuilderType<Patch, ValuesOnPatch<Patch>>::type;

    template <class Patch>
    using SplineOnPatch = DField<
            typename BuilderOnPatch<Patch>::batched_spline_domain_type<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    template <class Patch>
    using Derivs1OnPatch = DConstField<
            typename BuilderOnPatch<Patch>::batched_derivs_domain_type1<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    template <class Patch>
    using Derivs2OnPatch = DConstField<
            typename BuilderOnPatch<Patch>::batched_derivs_domain_type2<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    template <class Patch>
    using Derivs12OnPatch = DConstField<
            typename BuilderOnPatch<Patch>::batched_derivs_domain_type<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    using MultipatchSplineCoeffs = MultipatchField<SplineOnPatch, Patches...>;

    // For PERIODIC or GREVILLE boundary conditions
    using MultipatchValues = MultipatchField<ValuesOnPatch, Patches...>;

    using MultipatchDerivs1 = MultipatchField<Derivs1OnPatch, Patches...>;

    using MultipatchDerivs2 = MultipatchField<Derivs2OnPatch, Patches...>;

    using MultipatchDerivs12 = MultipatchField<Derivs12OnPatch, Patches...>;

    using BuilderTuple = std::tuple<BuilderOnPatch<Patches> const&...>;


    BuilderTuple const m_builders;

private:
    template <class Patch, template <typename P> typename DerivTypeOnPatch>
    std::optional<DerivTypeOnPatch<Patch>> get_deriv_value(
            std::optional<MultipatchField<DerivTypeOnPatch, Patches...>> derivs) const
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
            std::optional<MultipatchDerivs12> mixed_derivs_min1_min2 = std::nullopt,
            std::optional<MultipatchDerivs12> mixed_derivs_max1_min2 = std::nullopt,
            std::optional<MultipatchDerivs12> mixed_derivs_min1_max2 = std::nullopt,
            std::optional<MultipatchDerivs12> mixed_derivs_max1_max2 = std::nullopt) const
    {
        ((std::get<BuilderOnPatch<Patches> const&>(m_builders)(
                 splines.template get<Patches>(),
                 get_const_field(values.template get<Patches>()),
                 get_deriv_value<Patches, Derivs1OnPatch>(derivs_min1),
                 get_deriv_value<Patches, Derivs1OnPatch>(derivs_max1),
                 get_deriv_value<Patches, Derivs2OnPatch>(derivs_min2),
                 get_deriv_value<Patches, Derivs2OnPatch>(derivs_max2),
                 get_deriv_value<Patches, Derivs12OnPatch>(mixed_derivs_min1_min2),
                 get_deriv_value<Patches, Derivs12OnPatch>(mixed_derivs_max1_min2),
                 get_deriv_value<Patches, Derivs12OnPatch>(mixed_derivs_min1_max2),
                 get_deriv_value<Patches, Derivs12OnPatch>(mixed_derivs_max1_max2))),
         ...);
    };
};
```


