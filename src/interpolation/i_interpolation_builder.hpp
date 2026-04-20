// SPDX-License-Identifier: MIT
#pragma once

#include <optional>
#include <type_traits>

#include "ddc_aliases.hpp"

/**
 * @brief A traits struct for accessing type aliases of an interpolation builder.
 *
 * The primary template delegates to the builder's own type aliases, so any class
 * that defines them directly (e.g. IdentityInterpolationBuilder)
 * satisfies the InterpolationBuilder concept without specialisation.
 *
 * Specialise this struct to adapt external builders whose alias names differ from
 * the convention (e.g. ddc::SplineBuilder).
 *
 * Defines:
 *   Type aliases:
 *   - data_type
 *   - interpolation_idx_range_type
 *   - coeff_idx_range_type
 *   Static functions:
 *   - rank()
 *   Type calculators:
 *   - batched_basis_idx_range_type
 *   - batched_derivs_idx_range_type  (1D builders only; not required by the concept)
 *
 * @tparam Builder The interpolation builder type.
 */
template <class Builder>
struct InterpolationBuilderTraits
{
    /// @brief The data type that the data is saved on.
    using data_type = typename Builder::data_type;

    /// @brief The ND index range for the interpolation mesh.
    using interpolation_idx_range_type = typename Builder::interpolation_idx_range_type;

    /// @brief The index range for the interpolation coefficients.
    using coeff_idx_range_type = typename Builder::coeff_idx_range_type;

    /// @brief The number of interpolation dimensions.
    static constexpr std::size_t rank()
    {
        return interpolation_idx_range_type::rank();
    }

    /// @brief Batched domain with the interpolation grid(s) replaced by the basis grid(s).
    template <class IdxRangeBatchedInterpolation>
    using batched_basis_idx_range_type =
            typename Builder::template batched_basis_idx_range_type<IdxRangeBatchedInterpolation>;

    /// @brief Batched domain with the interpolation grid replaced by the deriv type.
    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type =
            typename Builder::template batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>;
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
    /// @brief The data type that the data is saved on.
    using data_type = double;

    /// @brief The discrete grid on which interpolation values are given.
    using interpolation_grid_type = typename Builder::interpolation_discrete_dimension_type;

    /// @brief The 1D index range for the interpolation mesh.
    using interpolation_idx_range_type = typename Builder::interpolation_domain_type;

    /// @brief The index range for the interpolation coefficients.
    using coeff_idx_range_type = IdxRange<typename Builder::bsplines_type>;

    /// @brief The number of interpolation dimensions (always 1 for SplineBuilder).
    static constexpr std::size_t rank()
    {
        return 1;
    }

    /// @brief The discrete dimension for the B-spline coefficients.
    using basis_domain_type = typename Builder::bsplines_type;

    /// @brief Batched domain with InterpolationDDim replaced by BSplines.
    template <class IdxRangeBatchedInterpolation>
    using batched_basis_idx_range_type =
            typename Builder::template batched_spline_domain_type<IdxRangeBatchedInterpolation>;

    /// @brief Batched domain with InterpolationDDim replaced by deriv_type.
    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type =
            typename Builder::template batched_derivs_domain_type<IdxRangeBatchedInterpolation>;
};

/**
 * @brief Get the batched basis index range for a builder.
 *
 * Dispatches to batched_basis_idx_range if available (e.g. IdentityInterpolationBuilder),
 * otherwise falls back to batched_spline_domain (e.g. ddc::SplineBuilder).
 *
 * @param builder The interpolation builder.
 * @param batched_interpolation_domain The batched interpolation domain.
 * @return The batched basis index range.
 */
template <class Builder, class IdxRangeBatchedInterpolation>
auto batched_basis_idx_range(
        Builder const& builder,
        IdxRangeBatchedInterpolation const& batched_interpolation_domain)
{
    if constexpr (requires { builder.batched_spline_domain(batched_interpolation_domain); }) {
        return builder.batched_spline_domain(batched_interpolation_domain);
    } else {
        return builder.batched_basis_idx_range(batched_interpolation_domain);
    }
}

namespace concepts {

/**
 * @brief A concept describing an ND interpolation builder.
 *
 * An interpolation builder is a callable that takes function values on an interpolation
 * mesh and computes the coefficients of an interpolation representation (e.g. spline or
 * Lagrange) of that function. The builder may operate over one or more interpolation
 * dimensions simultaneously (N ≥ 1).
 *
 * Type information is accessed through InterpolationBuilderTraits<Builder>, which has a
 * primary template delegating to Builder's own aliases and can be specialised for
 * external builders (e.g. ddc::SplineBuilder).
 *
 * The callable requirement is verified using interpolation_idx_range_type as the
 * representative non-batched domain:
 *  - b(coeffs, vals) — build interpolation coefficients from function values.
 */
template <class Builder>
concept InterpolationBuilder = requires
{
    typename Builder::exec_space;
    typename Builder::memory_space;
    typename InterpolationBuilderTraits<Builder>::data_type;
    typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type;
}
&&requires(
        Builder const& b,
        Field<typename InterpolationBuilderTraits<Builder>::data_type,
              typename InterpolationBuilderTraits<Builder>::template batched_basis_idx_range_type<
                      typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type>,
              typename Builder::memory_space> coeffs,
        ConstField<
                typename InterpolationBuilderTraits<Builder>::data_type,
                typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type,
                typename Builder::memory_space> vals)
{
    {b(coeffs, vals)};
};

/**
 * @brief A concept describing a 1D interpolation builder.
 *
 * Refines InterpolationBuilder with the additional requirements that:
 *   - The builder operates over exactly one interpolation dimension (rank() == 1).
 *   - InterpolationBuilderTraits<Builder>::basis_domain_type is defined, i.e. the
 *     builder exposes a discrete dimension for its basis coefficients.
 */
template <class Builder>
concept InterpolationBuilder1D
        = InterpolationBuilder<Builder> &&(InterpolationBuilderTraits<Builder>::rank() == 1)
          && requires(
                  Builder const& b,
                  typename InterpolationBuilderTraits<Builder>::interpolation_idx_range_type domain,
                  Field<typename InterpolationBuilderTraits<Builder>::data_type,
                        typename InterpolationBuilderTraits<Builder>::
                                template batched_basis_idx_range_type<
                                        typename InterpolationBuilderTraits<
                                                Builder>::interpolation_idx_range_type>,
                        typename Builder::memory_space> coeffs,
                  ConstField<
                          typename InterpolationBuilderTraits<Builder>::data_type,
                          typename InterpolationBuilderTraits<
                                  Builder>::interpolation_idx_range_type,
                          typename Builder::memory_space> vals,
                  std::optional<ConstField<
                          typename InterpolationBuilderTraits<Builder>::data_type,
                          typename InterpolationBuilderTraits<Builder>::
                                  template batched_derivs_idx_range_type<
                                          typename InterpolationBuilderTraits<
                                                  Builder>::interpolation_idx_range_type>,
                          typename Builder::memory_space>> derivs)
{
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
};

} // namespace concepts

/// @brief The discrete grid on which interpolation values are given (1D builders only).
template <concepts::InterpolationBuilder1D BuilderType>
using interpolation_grid_t = ddc::type_seq_element_t<
        0,
        ddc::to_type_seq_t<
                typename InterpolationBuilderTraits<BuilderType>::interpolation_idx_range_type>>;
