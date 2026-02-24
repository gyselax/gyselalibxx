// SPDX-License-Identifier: MIT
#pragma once

#include <optional>
#include <type_traits>

#include "ddc_aliases.hpp"

/**
 * @brief A traits struct for accessing type aliases of an interpolation builder.
 *
 * The primary template delegates to the builder's own type aliases, so any class
 * that defines them directly (e.g. IdentityInterpolationBuilder, SplineBuilder1D)
 * satisfies the InterpolationBuilder concept without specialisation.
 *
 * Specialise this struct to adapt external builders whose alias names differ from
 * the convention (e.g. ddc::SplineBuilder).
 *
 * Required type aliases in the primary template (or a specialisation):
 *   - exec_space                   : Kokkos execution space
 *   - memory_space                 : Kokkos memory space
 *   - continuous_dimension_type    : the continuous dimension being interpolated
 *   - interpolation_grid_type      : the discrete grid on which values are given
 *   - interpolation_idx_range_type : 1D index range for the interpolation mesh
 *   - basis_domain_type            : discrete dimension for the interpolation coefficients
 *   - deriv_type                   : ddc::Deriv tag for the boundary derivative dimension
 *   - batched_basis_idx_range_type<D>  : D with interpolation_grid_type replaced by basis_domain_type
 *   - batched_derivs_idx_range_type<D> : D with interpolation_grid_type replaced by deriv_type
 *
 * @tparam Builder The interpolation builder type.
 */
template <class Builder>
struct InterpolationBuilderTraits
{
    /// @brief The discrete grid on which interpolation values are given.
    using interpolation_grid_type = typename Builder::interpolation_grid_type;

    /// @brief The 1D index range for the interpolation mesh.
    using interpolation_idx_range_type = typename Builder::interpolation_idx_range_type;

    /// @brief The discrete dimension for the interpolation coefficients.
    using basis_domain_type = typename Builder::basis_domain_type;

    /// @brief Batched domain with interpolation_grid_type replaced by basis_domain_type.
    template <class BatchedInterpolationIdxRange>
    using batched_basis_idx_range_type =
            typename Builder::template batched_basis_idx_range_type<BatchedInterpolationIdxRange>;

    /// @brief Batched domain with interpolation_grid_type replaced by deriv_type.
    template <class BatchedInterpolationIdxRange>
    using batched_derivs_idx_range_type =
            typename Builder::template batched_derivs_idx_range_type<BatchedInterpolationIdxRange>;
};

/**
 * @brief Specialisation of InterpolationBuilderTraits for ddc::SplineBuilder.
 *
 * ddc::SplineBuilder uses different alias names from the InterpolationBuilder
 * convention. This specialisation provides the mapping so that ddc::SplineBuilder
 * can be used directly as an InterpolationBuilder without wrapping it in
 * SplineBuilder1D.
 *
 * Mapping:
 *   interpolation_discrete_dimension_type -> interpolation_grid_type
 *   interpolation_domain_type             -> interpolation_idx_range_type
 *   bsplines_type                         -> basis_domain_type
 *   batched_spline_domain_type<D>         -> batched_basis_idx_range_type<D>
 *   batched_derivs_domain_type<D>         -> batched_derivs_idx_range_type<D>
 */
template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class InterpolationDDim,
        ddc::BoundCond BcLower,
        ddc::BoundCond BcUpper,
        ddc::SplineSolver Solver>
struct InterpolationBuilderTraits<ddc::SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        InterpolationDDim,
        BcLower,
        BcUpper,
        Solver>>
{
private:
    using Builder = ddc::SplineBuilder<
            ExecSpace,
            MemorySpace,
            BSplines,
            InterpolationDDim,
            BcLower,
            BcUpper,
            Solver>;

public:
    /// @brief The discrete grid on which interpolation values are given.
    using interpolation_grid_type = typename Builder::interpolation_discrete_dimension_type;

    /// @brief The 1D index range for the interpolation mesh.
    using interpolation_idx_range_type = typename Builder::interpolation_domain_type;

    /// @brief The discrete dimension for the B-spline coefficients.
    using basis_domain_type = typename Builder::bsplines_type;

    /// @brief Batched domain with InterpolationDDim replaced by BSplines.
    template <class BatchedInterpolationIdxRange>
    using batched_basis_idx_range_type =
            typename Builder::template batched_spline_domain_type<BatchedInterpolationIdxRange>;

    /// @brief Batched domain with InterpolationDDim replaced by deriv_type.
    template <class BatchedInterpolationIdxRange>
    using batched_derivs_idx_range_type =
            typename Builder::template batched_derivs_domain_type<BatchedInterpolationIdxRange>;
};

/**
 * @brief Get the batched basis index range for a builder.
 *
 * Dispatches to @c batched_basis_idx_range if available (e.g. IdentityInterpolationBuilder),
 * otherwise falls back to @c batched_spline_domain (e.g. ddc::SplineBuilder).
 *
 * @param builder The interpolation builder.
 * @param batched_interpolation_domain The batched interpolation domain.
 * @return The batched basis index range.
 */
template <class Builder, class BatchedInterpolationIdxRange>
auto batched_basis_idx_range(
        Builder const& builder,
        BatchedInterpolationIdxRange const& batched_interpolation_domain)
{
    if constexpr (requires { builder.batched_spline_domain(batched_interpolation_domain); }) {
        return builder.batched_spline_domain(batched_interpolation_domain);
    } else {
        return builder.batched_basis_idx_range(batched_interpolation_domain);
    }
}

namespace concepts {

/**
 * @brief A concept describing an interpolation builder.
 *
 * An interpolation builder is a callable that takes function values on an interpolation
 * mesh and computes the coefficients of an interpolation representation (e.g. spline or
 * Lagrange) of that function.
 *
 * Type information is accessed through InterpolationBuilderTraits<Builder>, which has a
 * primary template delegating to Builder's own aliases and can be specialised for
 * external builders (e.g. ddc::SplineBuilder).
 *
 * A type satisfies this concept if InterpolationBuilderTraits<Builder> exposes:
 *
 * Non-template type aliases:
 *   - exec_space                   : Kokkos execution space
 *   - memory_space                 : Kokkos memory space
 *   - continuous_dimension_type    : the continuous dimension being interpolated
 *   - interpolation_grid_type      : the discrete grid on which values are given
 *   - interpolation_idx_range_type : 1D index range for the interpolation mesh
 *   - basis_domain_type            : discrete dimension for the interpolation coefficients
 *   - deriv_type                   : ddc::Deriv tag for the boundary derivative dimension
 *
 * Template type aliases (parameterised by any batched interpolation domain D):
 *   - batched_basis_idx_range_type<D>   : D with interpolation_grid_type replaced by basis_domain_type
 *   - batched_derivs_idx_range_type<D>  : D with interpolation_grid_type replaced by deriv_type
 *
 * Member functions (for a const builder b, a batched domain dom of type D, and matching fields):
 *   - b(coeffs, vals)                          : compute coefficients from values
 *   - b(coeffs, vals, derivs_min, derivs_max)  : compute coefficients with boundary derivatives
 *   - b.batched_derivs_xmin_domain(dom) -> batched_derivs_idx_range_type<D>
 *   - b.batched_derivs_xmax_domain(dom) -> batched_derivs_idx_range_type<D>
 *
 * The template alias and method requirements are verified using interpolation_idx_range_type
 * as a representative domain.
 */
template <class Builder>
concept InterpolationBuilder = requires
{
    typename Builder::exec_space;
    typename Builder::memory_space;
    typename Builder::continuous_dimension_type;
    typename InterpolationBuilderTraits<Builder>::interpolation_grid_type;
    typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type;
    typename InterpolationBuilderTraits<Builder>::basis_domain_type;
    typename Builder::deriv_type;
    // Verify the template aliases can be instantiated with a concrete domain
    typename InterpolationBuilderTraits<Builder>::template batched_basis_idx_range_type<
            typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type>;
    typename InterpolationBuilderTraits<Builder>::template batched_derivs_idx_range_type<
            typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type>;
}
&&requires(
        Builder const& b,
        typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type domain,
        Field<double,
              typename InterpolationBuilderTraits<Builder>::template batched_basis_idx_range_type<
                      typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type>,
              typename Builder::memory_space> coeffs,
        ConstField<
                double,
                typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type,
                typename Builder::memory_space> vals,
        std::optional<ConstField<
                double,
                typename InterpolationBuilderTraits<Builder>::
                        template batched_derivs_idx_range_type<typename InterpolationBuilderTraits<
                                Builder>::interpolation_idx_range_type>,
                typename Builder::memory_space>> derivs)
{
    {b(coeffs, vals)};
    {b(coeffs, vals, derivs, derivs)};
    {
        b.batched_derivs_xmin_domain(domain)
        } -> std::same_as<
                typename InterpolationBuilderTraits<Builder>::
                        template batched_derivs_idx_range_type<typename InterpolationBuilderTraits<
                                Builder>::interpolation_idx_range_type>>;
    {
        b.batched_derivs_xmax_domain(domain)
        } -> std::same_as<
                typename InterpolationBuilderTraits<Builder>::
                        template batched_derivs_idx_range_type<typename InterpolationBuilderTraits<
                                Builder>::interpolation_idx_range_type>>;
    {
        batched_basis_idx_range(b, domain)
        } -> std::same_as<
                typename InterpolationBuilderTraits<Builder>::template batched_basis_idx_range_type<
                        typename InterpolationBuilderTraits<
                                Builder>::interpolation_idx_range_type>>;
};

} // namespace concepts
