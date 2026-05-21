

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
KOKKOS_FUNCTION ddc::Real norm_inf(Coord<Tags...> coord)
{
    ddc::Real result = 0.0;
    ((result = Kokkos::max(result, Kokkos::fabs(coord.template get<Tags>()))), ...);
    return result;
}

template <class ElementType, class... Tags>
KOKKOS_FUNCTION ElementType norm_inf(Vector<ElementType, Tags...> vec)
{
    using index_set = typename Vector<ElementType, Tags...>::template vector_index_set_t<0>;
    static_assert(
            std::is_same_v<index_set, vector_index_set_dual_t<index_set>>,
            "Mapping is needed to calculate norm_inf on a non-orthonormal coordinate system");
    ElementType result = 0.0;
    ((result = Kokkos::max(result, Kokkos::fabs(ddcHelper::get<Tags>(vec)))), ...);
    return result;
}

template <std::floating_point T>
KOKKOS_INLINE_FUNCTION T norm_inf(T const coord)
{
    return Kokkos::fabs(coord);
}

namespace detail {

// General implementation of the infinity norm. This function is in a namespace to avoid code duplication
// without creating a function so general that it also captures multipatch types.
template <class DataType, class ExecSpace, class FuncType>
DataType norm_inf(ExecSpace exec_space, FuncType function)
{
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, typename FuncType::memory_space>::accessible);
    using IdxRangeFunc = typename FuncType::discrete_domain_type;
    using IdxFunc = typename IdxRangeFunc::discrete_element_type;
    IdxRangeFunc idx_range = get_idx_range(function);
    const std::source_location location = std::source_location::current();
    return ddc::parallel_transform_reduce(
            location.function_name(),
            exec_space,
            idx_range,
            0.,
            ddc::reducer::max<DataType>(),
            KOKKOS_LAMBDA(IdxFunc const idx) { return ::norm_inf(function(idx)); });
}

// General implementation of the infinity norm of an error. This function is in a namespace to avoid code duplication
// without creating a function so general that it also captures multipatch types.
template <class DataType, class ExecSpace, class FuncType, class ExactFuncType>
DataType error_norm_inf(ExecSpace exec_space, FuncType function, ExactFuncType exact_function)
{
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, typename FuncType::memory_space>::accessible);
    using IdxRangeFunc = typename FuncType::discrete_domain_type;
    using IdxFunc = typename IdxRangeFunc::discrete_element_type;
    IdxRangeFunc idx_range = get_idx_range(function);
    const std::source_location location = std::source_location::current();
    return ddc::parallel_transform_reduce(
            location.function_name(),
            exec_space,
            idx_range,
            DataType(0.),
            ddc::reducer::max<DataType>(),
            KOKKOS_LAMBDA(IdxFunc const idx) {
                return ::norm_inf(function(idx) - exact_function(idx));
            });
}

}; // namespace detail

template <class ExecSpace, class ElementType, class IdxRange>
inline ElementType norm_inf(
        ExecSpace exec_space,
        ConstField<ElementType, IdxRange, typename ExecSpace::memory_space> function)
{
    return detail::norm_inf<ElementType>(exec_space, function);
}

template <class ExecSpace, class IdxRange, class... CoordDims>
inline ddc::Real norm_inf(
        ExecSpace exec_space,
        ConstField<Coord<CoordDims...>, IdxRange, typename ExecSpace::memory_space> function)
{
    return detail::norm_inf<ddc::Real>(exec_space, function);
}

template <class ExecSpace, class ElementType, class IdxRange, class VectorIndexSetType>
inline ElementType norm_inf(
        ExecSpace exec_space,
        VectorConstField<
                ElementType,
                IdxRange,
                VectorIndexSetType,
                typename ExecSpace::memory_space> function)
{
    return detail::norm_inf<ElementType>(exec_space, function);
}

template <class ExecSpace, class ElementType, class IdxRange>
inline ElementType error_norm_inf(
        ExecSpace exec_space,
        ConstField<ElementType, IdxRange, typename ExecSpace::memory_space> function,
        ConstField<ElementType, IdxRange, typename ExecSpace::memory_space> exact_function)
{
    return detail::error_norm_inf<ElementType>(exec_space, function, exact_function);
}

template <class ExecSpace, class IdxRange, class... CoordDims>
inline ddc::Real error_norm_inf(
        ExecSpace exec_space,
        ConstField<Coord<CoordDims...>, IdxRange, typename ExecSpace::memory_space> function,
        ConstField<Coord<CoordDims...>, IdxRange, typename ExecSpace::memory_space> exact_function)
{
    return detail::error_norm_inf<ddc::Real>(exec_space, function, exact_function);
}

namespace concepts {
template <class F, class ElementType, class IdxType>
concept CompatibleFunc = requires(F const& f, IdxType idx)
{
    {
        f(idx)
        } -> std::same_as<ElementType>;
};
} // namespace concepts

template <class ExecSpace, class ElementType, class IdxRange, class ExactFunc>
inline ElementType error_norm_inf(
        ExecSpace exec_space,
        ConstField<ElementType, IdxRange, typename ExecSpace::memory_space> function,
        ExactFunc exact_function)
{
    static_assert(concepts::CompatibleFunc<
                  ExactFunc,
                  ElementType,
                  typename IdxRange::discrete_element_type>);
    return detail::error_norm_inf<ElementType>(exec_space, function, exact_function);
}

template <class ExecSpace, class ElementType, class IdxRange, class VectorIndexSetType>
inline ElementType error_norm_inf(
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
    return detail::error_norm_inf<ElementType>(exec_space, function, exact_function);
}

template <class IdxRangeQuad, class ExecSpace, class DataType>
DataType norm_L1(
        ExecSpace exec_space,
        Quadrature<IdxRangeQuad, IdxRangeQuad, DataType, typename ExecSpace::memory_space>
                quadrature,
        Field<DataType, IdxRangeQuad, typename ExecSpace::memory_space> function)
{
    using IdxQuad = typename IdxRangeQuad::discrete_element_type;
    return quadrature(
            exec_space,
            KOKKOS_LAMBDA(IdxQuad const idx) { return Kokkos::fabs(function(idx)); });
}

template <class IdxRangeQuad, class ExecSpace, class DataType>
DataType error_norm_L1(
        ExecSpace exec_space,
        Quadrature<IdxRangeQuad, IdxRangeQuad, DataType, typename ExecSpace::memory_space>
                quadrature,
        Field<DataType, IdxRangeQuad, typename ExecSpace::memory_space> function,
        Field<DataType, IdxRangeQuad, typename ExecSpace::memory_space> exact_function)
{
    using IdxQuad = typename IdxRangeQuad::discrete_element_type;
    return quadrature(
            exec_space,
            KOKKOS_LAMBDA(IdxQuad const idx) {
                return Kokkos::fabs(function(idx) - exact_function(idx));
            });
}


template <class IdxRangeQuad, class ExecSpace, class DataType>
DataType norm_L2(
        ExecSpace exec_space,
        Quadrature<IdxRangeQuad, IdxRangeQuad, DataType, typename ExecSpace::memory_space>
                quadrature,
        Field<DataType, IdxRangeQuad, typename ExecSpace::memory_space> function)
{
    using IdxQuad = typename IdxRangeQuad::discrete_element_type;
    return std::sqrt(quadrature(
            exec_space,
            KOKKOS_LAMBDA(IdxQuad const idx) { return function(idx) * function(idx); }));
}

template <class IdxRangeQuad, class ExecSpace, class DataType>
DataType error_norm_L2(
        ExecSpace exec_space,
        Quadrature<IdxRangeQuad, IdxRangeQuad, DataType, typename ExecSpace::memory_space>
                quadrature,
        Field<DataType, IdxRangeQuad, typename ExecSpace::memory_space> function,
        Field<DataType, IdxRangeQuad, typename ExecSpace::memory_space> exact_function)
{
    using IdxQuad = typename IdxRangeQuad::discrete_element_type;
    return std::sqrt(quadrature(
            exec_space,
            KOKKOS_LAMBDA(IdxQuad const idx) {
                DataType err = function(idx) - exact_function(idx);
                return err * err;
            }));
}
```


