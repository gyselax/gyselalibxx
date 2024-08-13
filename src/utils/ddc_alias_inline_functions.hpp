// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>

#include <ddc/ddc.hpp>

#include <sll/polar_spline.hpp>

#include "derivative_field_common.hpp"
#include "vector_field_common.hpp"

/**
 * A function to get the range of valid indices that can be used to index this field.
 *
 * @param[in] field The field whose indices are of interest.
 *
 * @returns The index range.
 */
template <class... QueryGrids, class FieldType>
KOKKOS_INLINE_FUNCTION auto get_idx_range(FieldType const& field) noexcept
{
    static_assert(
            ddc::is_chunk_v<FieldType> || is_field_v<FieldType>,
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

/**
 * A function to get the range of valid indices that can be used to index a set of b-splines that are compatible with this spline builder.
 *
 * @param[in] builder The spline builder.
 *
 * @returns The index range.
 */
template <class SplineBuilder>
inline auto get_spline_idx_range(SplineBuilder const& builder) noexcept
{
    return builder.spline_domain();
}

/**
 * A helper function to get a modifiable field from a FieldMem without allocating additional memory.
 *
 * @param[in] field The field memory object.
 * @returns The modifiable field.
 */
template <class FieldType>
inline auto get_field(FieldType&& field)
{
    using Type = std::remove_cv_t<std::remove_reference_t<FieldType>>;
    static_assert(
            (ddc::is_chunk_v<Type>) || (is_field_v<Type>) || (is_deriv_field_v<Type>)
                    || (is_polar_spline_v<Type>),
            "Not a Field or FieldMem (ddc::Chunk or ddc::ChunkSpan) type");
    return field.span_view();
}

/**
 * A helper function to get a constant field from a FieldMem without allocating additional memory.
 *
 * @param[in] field The field memory object.
 * @returns The constant field.
 */
template <class FieldType>
inline auto get_const_field(FieldType&& field)
{
    using Type = std::remove_cv_t<std::remove_reference_t<FieldType>>;
    static_assert(
            (ddc::is_chunk_v<Type>) || (is_field_v<Type>) || (is_deriv_field_v<Type>)
                    || (is_polar_spline_v<Type>),
            "Not a Field or FieldMem (ddc::Chunk or ddc::ChunkSpan) type");
    return field.span_cview();
}
