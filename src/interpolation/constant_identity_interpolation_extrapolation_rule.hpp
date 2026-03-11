// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_aliases.hpp"

/**
 * @brief A constant extrapolation rule for identity-based (Lagrange) interpolation.
 *
 * When an evaluation coordinate falls outside the interpolation domain, this rule
 * returns the interpolation coefficient at a fixed boundary index (@c m_coeff_idx)
 * rather than extrapolating the polynomial. This is the Lagrange analog of
 * ddc::ConstantExtrapolationRule used with spline evaluators.
 *
 * It is intended to be constructed via get_extrapolation<CONSTANT, CoeffGrid, Basis>(),
 * which supplies the appropriate boundary index automatically.
 *
 * @tparam CoeffGrid The discrete grid type of the interpolation coefficients.
 * @tparam DataType  The floating-point type of the function values (default: double).
 */
template <class CoeffGrid, class DataType = double>
class ConstantIdentityInterpolationExtrapolationRule
{
private:
    Idx<CoeffGrid> m_coeff_idx;

public:
    /**
     * @brief Construct the rule, clamping to the coefficient at @c coeff_idx.
     *
     * @param coeff_idx The index of the boundary coefficient to return for all
     *                  out-of-domain evaluation coordinates.
     */
    explicit ConstantIdentityInterpolationExtrapolationRule(Idx<CoeffGrid> coeff_idx)
        : m_coeff_idx(coeff_idx)
    {
    }

    /**
     * @brief Return the boundary interpolation coefficient for a coordinate outside the domain.
     *
     * @param[in] pos         The coordinate where the function is to be evaluated (unused).
     * @param[in] interp_coef The full field of interpolation coefficients.
     *
     * @return The value of @c interp_coef at the fixed boundary index.
     */
    template <class CoordType, class Layout, class MemorySpace>
    KOKKOS_FUNCTION double operator()(
            [[maybe_unused]] CoordType pos,
            ConstField<DataType, IdxRange<CoeffGrid>, Layout, MemorySpace> const interp_coef) const
    {
        return interp_coef(m_coeff_idx);
    }
};
