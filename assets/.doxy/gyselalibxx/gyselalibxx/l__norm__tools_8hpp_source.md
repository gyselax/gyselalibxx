

# File l\_norm\_tools.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**l\_norm\_tools.hpp**](l__norm__tools_8hpp.md)

[Go to the documentation of this file](l__norm__tools_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "quadrature.hpp"
#include "vector_field.hpp"
#include "vector_index_tools.hpp"

template <class... Tags>
KOKKOS_FUNCTION double norm_inf(ddc::Coordinate<Tags...> coord)
{
    double result = 0.0;
    ((result = Kokkos::max(result, Kokkos::fabs(coord.template get<Tags>()))), ...);
    return result;
}

template <class... Tags>
KOKKOS_FUNCTION double norm_inf(DVector<Tags...> vec)
{
    using index_set = typename DVector<Tags...>::template vector_index_set_t<0>;
    static_assert(
            std::is_same_v<index_set, vector_index_set_dual_t<index_set>>,
            "Mapping is needed to calculate norm_inf on a non-orthonormal coordinate system");
    double result = 0.0;
    ((result = Kokkos::max(result, Kokkos::fabs(ddcHelper::get<Tags>(vec)))), ...);
    return result;
}

KOKKOS_INLINE_FUNCTION double norm_inf(double const coord)
{
    return Kokkos::fabs(coord);
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

template <class ExecSpace, class ElementType, class IdxRange>
inline double norm_inf(
        ExecSpace exec_space,
        ConstField<ElementType, IdxRange, typename ExecSpace::memory_space> function)
{
    return detail::norm_inf(exec_space, function);
}

template <class ExecSpace, class ElementType, class IdxRange, class VectorIndexSetType>
inline double norm_inf(
        ExecSpace exec_space,
        VectorConstField<
                ElementType,
                IdxRange,
                VectorIndexSetType,
                typename ExecSpace::memory_space> function)
{
    return detail::norm_inf(exec_space, function);
}

template <class ExecSpace, class ElementType, class IdxRange>
inline double error_norm_inf(
        ExecSpace exec_space,
        ConstField<ElementType, IdxRange, typename ExecSpace::memory_space> function,
        ConstField<ElementType, IdxRange, typename ExecSpace::memory_space> exact_function)
{
    return detail::error_norm_inf(exec_space, function, exact_function);
}

template <class ExecSpace, class ElementType, class IdxRange, class VectorIndexSetType>
inline double error_norm_inf(
        ExecSpace exec_space,
        VectorConstField<
                ElementType,
                IdxRange,
                VectorIndexSetType,
                typename ExecSpace::memory_space> function,
        VectorConstField<
                ElementType,
                IdxRange,
                VectorIndexSetType,
                typename ExecSpace::memory_space> exact_function)
{
    return detail::error_norm_inf(exec_space, function, exact_function);
}

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
```


