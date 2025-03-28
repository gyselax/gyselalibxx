

# File ddc\_alias\_inline\_functions.hpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**ddc\_alias\_inline\_functions.hpp**](ddc__alias__inline__functions_8hpp.md)

[Go to the documentation of this file](ddc__alias__inline__functions_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>
#include <utility>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

namespace detail {

/*
 * A class using template magic to determine if a type has an idx_range function.
 */
template <typename Type>
class HasIdxRange
{
private:
    // Class for SFINAE deduction
    template <typename U>
    class check
    {
    };

    // Function that will be chosen if ClassType has a method called idx_range with 0 arguments
    template <typename ClassType>
    static char f(check<decltype(std::declval<ClassType>().idx_range())>*);

    // Function that will be chosen by default
    template <typename ClassType>
    static long f(...);

public:
    static constexpr bool value = (sizeof(f<Type>(0)) == sizeof(char));
};

/*
 * A class using template magic to determine if a type is a kind of Field defined in Gyselalibxx.
 * Such a field should define a function get_const_field.
 */
template <typename Type>
class IsGslxField
{
private:
    // Class for SFINAE deduction
    template <typename U>
    class check
    {
    };

    // Function that will be chosen if ClassType has a method called get_const_field with 0 arguments
    template <typename ClassType>
    static char f(check<decltype(std::declval<ClassType>().get_const_field())>*);

    // Function that will be chosen by default
    template <typename ClassType>
    static long f(...);

public:
    static constexpr bool value = (sizeof(f<Type>(0)) == sizeof(char));
};

} // namespace detail
//
template <class T>
inline constexpr bool enable_data_access_methods = false;

template <class T>
inline constexpr bool enable_mem_type = false;

template <class ElementType, class IdxRangeType, class LayoutType, class MemoryType>
inline constexpr bool
        enable_data_access_methods<Field<ElementType, IdxRangeType, LayoutType, MemoryType>> = true;

template <class ElementType, class IdxRangeType, class MemoryType>
inline constexpr bool
        enable_data_access_methods<FieldMem<ElementType, IdxRangeType, MemoryType>> = true;

template <class ElementType, class IdxRangeType, class MemoryType>
inline constexpr bool enable_mem_type<FieldMem<ElementType, IdxRangeType, MemoryType>> = true;

template <typename Type>
static constexpr bool has_idx_range_v = detail::HasIdxRange<Type>::value;

template <typename Type>
static constexpr bool is_gslx_field_v = detail::IsGslxField<Type>::value;

template <typename Type>
inline constexpr bool has_data_access_methods_v
        = enable_data_access_methods<std::remove_const_t<std::remove_reference_t<Type>>>;

template <class Type>
inline constexpr bool is_mem_type_v
        = enable_mem_type<std::remove_const_t<std::remove_reference_t<Type>>>;

template <class... QueryGrids, class FieldType>
auto get_idx_range(FieldType const& field) noexcept
{
    static_assert(
            ddc::is_chunk_v<FieldType> || has_idx_range_v<FieldType>,
            "Not a DDC field (ddc::ChunkSpan) type");
    if constexpr (ddc::is_chunk_v<FieldType>) {
        if constexpr (sizeof...(QueryGrids) == 0) {
            return field.domain();
        } else {
            return field.template domain<QueryGrids...>();
        }
    } else {
        if constexpr (sizeof...(QueryGrids) == 0) {
            return field.idx_range();
        } else {
            return field.template idx_range<QueryGrids...>();
        }
    }
}

template <class SplineBuilder>
inline auto get_spline_idx_range(SplineBuilder const& builder) noexcept
{
    return builder.spline_domain();
}

template <class FieldType, std::enable_if_t<!is_mem_type_v<FieldType>, bool> = true>
KOKKOS_INLINE_FUNCTION auto get_field(FieldType&& field)
{
    static_assert(
            has_data_access_methods_v<FieldType>,
            "Not a Field or FieldMem (ddc::Chunk or ddc::ChunkSpan) type");
    return field.span_view();
}

template <class FieldType, std::enable_if_t<is_mem_type_v<FieldType>, bool> = true>
inline auto get_field(FieldType&& field)
{
    static_assert(
            has_data_access_methods_v<FieldType>,
            "Not a Field or FieldMem (ddc::Chunk or ddc::ChunkSpan) type");
    return field.span_view();
}

template <class FieldType, std::enable_if_t<!is_mem_type_v<FieldType>, bool> = true>
KOKKOS_INLINE_FUNCTION auto get_const_field(FieldType&& field)
{
    static_assert(
            has_data_access_methods_v<FieldType>,
            "Not a Field or FieldMem (ddc::Chunk or ddc::ChunkSpan) type");
    return field.span_cview();
}

template <class FieldType, std::enable_if_t<is_mem_type_v<FieldType>, bool> = true>
inline auto get_const_field(FieldType&& field)
{
    static_assert(
            has_data_access_methods_v<FieldType>,
            "Not a Field or FieldMem (ddc::Chunk or ddc::ChunkSpan) type");
    return field.span_cview();
}
```


