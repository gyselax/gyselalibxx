// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>
#include <tuple>
#include <utility>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"

/**
 * @brief A class to call all the builders of all the patches once.
 *
 * We need to instantiate all the builders for all the patches in the main code.
 * We process the same way for the Fields containing the spline coefficients and the
 * values of the function on each patch. The Fields are stored in MultipatchField
 * objects. The builders are stored in this class.
 * This class is instantiated with all the builders.
 * The operator() allows to call all the builders stored in the member of this class in 
 * one single line.
 *
 * This function is useful to avoid calling all the builders individually, especially in
 * a multipatch geometry with several patches.
 *
 * @tparam ExecSpace The space (CPU/GPU) where the calculations are carried out.
 * @tparam MemorySpace The space (CPU/GPU) where the coefficients and values are stored.
 * @tparam BSplineOnPatch A type alias which provides the BSpline type along which the splines are built.
 * @tparam GridOnPatch A type alias which provides the Grid type along which the interpolation points of the splines are found.
 * @tparam BcLower The lower boundary condition.
 * @tparam BcUpper The upper boundary condition.
 * @tparam BcTransition The boundary condition used at the interface between 2 patches.
 * @tparam Connectivity A MultipatchConnectivity object describing the interfaces between patches.
 * @tparam Solver The SplineSolver giving the backend used to perform the spline approximation. See DDC for more details.
 * @tparam ValuesOnPatch A type alias which provides the field type which will be used to pass the values of the function
 *                      at the interpolation points. The index range of this field type should contain any batch dimensions.
 */
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

    /**
     * A small structure allowing the multiple grids to be unpacked from a field and repacked into
     * a SplineBuilder type.
     */
    template <class Patch>
    struct Build_BuilderType
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
                Solver>;
    };

    /// A type alias to get the builder type on a specific patch.
    template <class Patch>
    using BuilderOnPatch = typename Build_BuilderType<Patch>::type;

    /// A type alias to get the batched spline coefficients on a specific patch.
    template <class Patch>
    using SplineOnPatch = DField<
            typename BuilderOnPatch<Patch>::batched_spline_domain_type<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    /// A type alias to get the batched derivatives on a specific patch.
    template <class Patch>
    using DerivsOnPatch = DConstField<
            typename BuilderOnPatch<Patch>::batched_derivs_domain_type<
                    typename ValuesOnPatch<Patch>::discrete_domain_type>,
            MemorySpace>;

    /// The type of the batched spline coefficients.
    using MultipatchSplineCoeffs = MultipatchField<SplineOnPatch, Patches...>;

    // For PERIODIC or GREVILLE boundary conditions
    /// The type of the values at the batched interpolation points.
    using MultipatchValues = MultipatchField<ValuesOnPatch, Patches...>;

    using MultipatchDerivs = MultipatchField<DerivsOnPatch, Patches...>;

    /// The type of the internal storage of the SplineBuilders.
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
    /**
     * @brief Instantiate the MultipatchSplineBuilder from a std::tuple 
     * of all the builder on each patch. 
     * 
     * @warning The builders have to be sorted in the same order as the patches
     * in the tuple. 
     * 
     * @param builders Spline builders for each patch.
     */
    explicit MultipatchSplineBuilder(BuilderOnPatch<Patches> const&... builders)
        : m_builders(std::tie(builders...)) {};


    ~MultipatchSplineBuilder() = default;

    /**
     * @brief Build the spline representation of each given function.
     * 
     * @param splines MultipatchField of all the Fields pointing to the spline representations. 
     * @param values MultipatchField of all the Fields pointing to the function values. 
     * @param derivs_xmin MultipatchField of all the ConstFields describing the function derivatives at the lower bound.
     * @param derivs_xmax MultipatchField of all the ConstFields describing the function derivatives at the upper bound.
     */
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
