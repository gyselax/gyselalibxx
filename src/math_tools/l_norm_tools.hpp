// SPDX-License-Identifier: MIT
/**
 * @file l_norm_tools.hpp
 * File Describing useful mathematical functions to compute Lnorms
 */

#pragma once

#include <ddc/ddc.hpp>

#include "quadrature.hpp"
#include "vector_field.hpp"

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
 * @brief Compute the infinity norm of a vector on an
 * orthonormal coordinate system.
 *
 * For a given vector @f$ x @f$ , compute
 * @f$|Vert x |Vert_{\infty} = \sup_n |x_n| @f$.
 *
 * @param[in] vec The given vector.
 *
 * @return A double containing the value of the infinty norm.
 */
template <class... Tags>
KOKKOS_FUNCTION double norm_inf(DVector<Tags...> vec)
{
    using index_set = typename DVector<Tags...>::vector_index_set_t<0>;
    static_assert(
            std::is_same_v<index_set, vector_index_set_dual_t<index_set>>,
            "Mapping is needed to calculate norm_inf on a non-orthonormal coordinate system");
    double result = 0.0;
    ((result = Kokkos::max(result, Kokkos::fabs(ddcHelper::get<Tags>(vec)))), ...);
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

namespace detail {

// General implementation of the infinity norm. This function is in a namespace to avoid code duplication
// without creating a function so general that it also captures multipatch types.
template <class ExecSpace, class FuncType>
double norm_inf(ExecSpace exec_space, FuncType function)
{
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, typename FuncType::memory_space>::accessible);
    using IdxRangeFunc = typename FuncType::discrete_domain_type;
    using IdxFunc = typename IdxRangeFunc::discrete_element_type;
    IdxRangeFunc idx_range = get_idx_range(function);
    return ddc::parallel_transform_reduce(
            exec_space,
            idx_range,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxFunc const idx) { return ::norm_inf(function(idx)); });
}

// General implementation of the infinity norm of an error. This function is in a namespace to avoid code duplication
// without creating a function so general that it also captures multipatch types.
template <class ExecSpace, class FuncType>
double error_norm_inf(ExecSpace exec_space, FuncType function, FuncType exact_function)
{
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, typename FuncType::memory_space>::accessible);
    using IdxRangeFunc = typename FuncType::discrete_domain_type;
    using IdxFunc = typename IdxRangeFunc::discrete_element_type;
    IdxRangeFunc idx_range = get_idx_range(function);
    return ddc::parallel_transform_reduce(
            exec_space,
            idx_range,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxFunc const idx) {
                return ::norm_inf(function(idx) - exact_function(idx));
            });
}

}; // namespace detail

/**
 * @brief Compute the infinity norm for a Field.
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] function The function whose norm is calcuated.
 * @return A double containing the value of the infinty norm.
 */
template <class ExecSpace, class ElementType, class IdxRange>
inline double norm_inf(
        ExecSpace exec_space,
        ConstField<ElementType, IdxRange, typename ExecSpace::memory_space> function)
{
    return detail::norm_inf(exec_space, function);
}

/**
 * @brief Compute the infinity norm for a VectorField.
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] function The function whose norm is calcuated.
 * @return A double containing the value of the infinty norm.
 */
template <class ExecSpace, class ElementType, class IdxRange, class NDTag>
inline double norm_inf(
        ExecSpace exec_space,
        VectorConstField<ElementType, IdxRange, NDTag, typename ExecSpace::memory_space> function)
{
    return detail::norm_inf(exec_space, function);
}

/**
 * @brief Compute the infinity norm of the error between 2 Fields.
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] function The calculated function.
 * @param[in] exact_function The exact function with which the calculated function is compared.
 * @return A double containing the value of the infinty norm.
 */
template <class ExecSpace, class ElementType, class IdxRange>
inline double error_norm_inf(
        ExecSpace exec_space,
        ConstField<ElementType, IdxRange, typename ExecSpace::memory_space> function,
        ConstField<ElementType, IdxRange, typename ExecSpace::memory_space> exact_function)
{
    return detail::error_norm_inf(exec_space, function, exact_function);
}

/**
 * @brief Compute the infinity norm of the error between 2 VectorFields.
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] function The calculated function.
 * @param[in] exact_function The exact function with which the calculated function is compared.
 * @return A double containing the value of the infinty norm.
 */
template <class ExecSpace, class ElementType, class IdxRange, class NDTag>
inline double error_norm_inf(
        ExecSpace exec_space,
        VectorConstField<ElementType, IdxRange, NDTag, typename ExecSpace::memory_space> function,
        VectorConstField<ElementType, IdxRange, NDTag, typename ExecSpace::memory_space>
                exact_function)
{
    return detail::error_norm_inf(exec_space, function, exact_function);
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
        DField<IdxRangeQuad, typename ExecSpace::memory_space> function)
{
    using IdxQuad = typename IdxRangeQuad::discrete_element_type;
    return quadrature(
            exec_space,
            KOKKOS_LAMBDA(IdxQuad const idx) { return Kokkos::fabs(function(idx)); });
}

/**
 * @brief Compute the L1 norm of the error between 2 Fields.
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] quadrature The quadrature used to compute the integral.
 * @param[in] function The calculated function.
 * @param[in] exact_function The exact function with which the calculated function is compared.
 * @return A double containing the value of the infinty norm.
 */
template <class IdxRangeQuad, class ExecSpace>
double error_norm_L1(
        ExecSpace exec_space,
        Quadrature<IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory_space> quadrature,
        DField<IdxRangeQuad, typename ExecSpace::memory_space> function,
        DField<IdxRangeQuad, typename ExecSpace::memory_space> exact_function)
{
    using IdxQuad = typename IdxRangeQuad::discrete_element_type;
    return quadrature(
            exec_space,
            KOKKOS_LAMBDA(IdxQuad const idx) {
                return Kokkos::fabs(function(idx) - exact_function(idx));
            });
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

/**
 * @brief Compute the L2 norm of the error between 2 Fields.
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] quadrature The quadrature used to compute the integral.
 * @param[in] function The calculated function.
 * @param[in] exact_function The exact function with which the calculated function is compared.
 * @return A double containing the value of the infinty norm.
 */
template <class IdxRangeQuad, class ExecSpace>
double error_norm_L2(
        ExecSpace exec_space,
        Quadrature<IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory_space> quadrature,
        DField<IdxRangeQuad, typename ExecSpace::memory_space> function,
        DField<IdxRangeQuad, typename ExecSpace::memory_space> exact_function)
{
    using IdxQuad = typename IdxRangeQuad::discrete_element_type;
    return std::sqrt(quadrature(
            exec_space,
            KOKKOS_LAMBDA(IdxQuad const idx) {
                double err = function(idx) - exact_function(idx);
                return err * err;
            }));
}
