#pragma once
#include "ddc_aliases.hpp"

template <class Basis, class DataType = double>
class ConstantIdentityInterpolationExtrapolationRule
{
private:
    Idx<Basis> m_coeff_idx;

public:
    ConstantIdentityInterpolationExtrapolationRule(Idx<Basis> coeff_idx) : m_coeff_idx(coeff_idx) {}

    /**
     * @brief Get the value of the function on B-splines at a coordinate outside the domain.
     *
     * @param[in] pos The coordinate where we want to evaluate the function on B-splines.
     * @param[in] spline_coef The coefficients of the function on B-splines.
     *
     * @return A double with the value of the function on B-splines evaluated at the coordinate.
     */
    template <class CoordType, class Layout, class MemorySpace>
    KOKKOS_FUNCTION double operator()(
            [[maybe_unused]] CoordType pos,
            ConstField<DataType, IdxRange<Basis>, Layout, MemorySpace> const interp_coef) const
    {
        return interp_coef(m_coeff_idx);
    }
};
