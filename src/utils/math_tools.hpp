// SPDX-License-Identifier: MIT
/**
 * @file math_tools.hpp
 * File Describing useful mathematical functions.
 */

#pragma once

#include <ddc/ddc.hpp>

#include "quadrature.hpp"

/**
 * @brief Compute the infinity norm.
 *
 * For a given vector @f$ x @f$ , compute
 * @f$|Vert x |Vert_{\infty} = \sup_n |x_n| @f$.
 *
 * @param[in] coord
 *      The given vector.
 *
 * @return A double containing the value of the infinty norm.
 */
template <class... Tags>
KOKKOS_FUNCTION double norm_inf(ddc::Coordinate<Tags...> coord)
{
    double result = 0.0;
    ((result = Kokkos::max(result, Kokkos::fabs(coord.template get<Tags>()))), ...);
    return result;
}


/**
 * @brief Compute the infinity norm.
 *
 * In case of scalar, the infinity norm
 * returns the scalar.
 *
 * @param[in] coord
 *      The given double.
 *
 * @return A double containing the value of the infinty norm.
 */
KOKKOS_INLINE_FUNCTION double norm_inf(double const coord)
{
    return coord;
}

/**
 * @brief Compute L1 norm of a function with a given quadrature.
 *
 * @f$ \int_{\Omega} |f(X)|  dX @f$
 *
 * @param[in] exec_space
 *     The space on which the function is executed (CPU/GPU).
 * @param[in] quadrature
 *      The quadrature used to compute the integral.
 * @param[in] function
 *      A Field to the value of the function on the quadrature grid.
 *
 * @return A double containing the L1 norm of the function.
 */
template <class IdxRangeQuad, class ExecSpace>
double norm_L1(
        ExecSpace exec_space,
        Quadrature<IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory_space> quadrature,
        Field<double, IdxRangeQuad, typename ExecSpace::memory_space> function)
{
    using IdxQuad = typename IdxRangeQuad::discrete_element_type;
    return quadrature(
            exec_space,
            KOKKOS_LAMBDA(IdxQuad const idx) { return Kokkos::fabs(function(idx)); });
}



/**
 * @brief Compute L2 norm of a function with a given quadrature.
 *
 * @f$ \sqrt{\int_{\Omega} |f(X)|^2  dX} @f$
 *
 * @param[in] exec_space
 *     The space on which the function is executed (CPU/GPU).
 * @param[in] quadrature
 *      The quadrature used to compute the integral.
 * @param[in] function
 *      A Field to the value of the function on the quadrature grid.
 *
 * @return A double containing the L2 norm of the function.
 */
template <class IdxRangeQuad, class ExecSpace>
double norm_L2(
        ExecSpace exec_space,
        Quadrature<IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory_space> quadrature,
        DField<IdxRangeQuad, typename ExecSpace::memory_space> function)
{
    using IdxQuad = typename IdxRangeQuad::discrete_element_type;
    return std::sqrt(quadrature(
            exec_space,
            KOKKOS_LAMBDA(IdxQuad const idx) { return function(idx) * function(idx); }));
}
