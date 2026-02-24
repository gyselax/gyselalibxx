// SPDX-License-Identifier: MIT
#pragma once

#include <optional>
#include <type_traits>

#include "ddc_aliases.hpp"

namespace concepts {

/**
 * @brief A concept describing an interpolation builder.
 *
 * An interpolation builder is a callable that takes function values on an interpolation
 * mesh and computes the coefficients of an interpolation representation (e.g. spline or
 * Lagrange) of that function. A type satisfies this concept if it exposes:
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
    typename Builder::interpolation_grid_type;
    typename Builder::interpolation_idx_range_type;
    typename Builder::basis_domain_type;
    typename Builder::deriv_type;
    // Verify the template aliases can be instantiated with a concrete domain
    typename Builder::template batched_basis_idx_range_type<
            typename Builder::interpolation_idx_range_type>;
    typename Builder::template batched_derivs_idx_range_type<
            typename Builder::interpolation_idx_range_type>;
}
&&requires(
        Builder const& b,
        typename Builder::interpolation_idx_range_type domain,
        Field<double,
              typename Builder::template batched_basis_idx_range_type<
                      typename Builder::interpolation_idx_range_type>,
              typename Builder::memory_space> coeffs,
        ConstField<
                double,
                typename Builder::interpolation_idx_range_type,
                typename Builder::memory_space> vals,
        std::optional<ConstField<
                double,
                typename Builder::template batched_derivs_idx_range_type<
                        typename Builder::interpolation_idx_range_type>,
                typename Builder::memory_space>> derivs)
{
    {b(coeffs, vals)};
    {b(coeffs, vals, derivs, derivs)};
    {
        b.batched_derivs_xmin_domain(domain)
        } -> std::same_as<typename Builder::template batched_derivs_idx_range_type<
                typename Builder::interpolation_idx_range_type>>;
    {
        b.batched_derivs_xmax_domain(domain)
        } -> std::same_as<typename Builder::template batched_derivs_idx_range_type<
                typename Builder::interpolation_idx_range_type>>;
};

} // namespace concepts
