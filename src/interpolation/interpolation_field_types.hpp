// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_aliases.hpp"
#include "vector_field.hpp"

namespace detail {

/**
 * @brief The field types used to pass data to/from an interpolation builder or evaluator
 * acting on a single index range.
 *
 * This struct is used as a base class of InterpolationBuilderTraits and
 * InterpolationEvaluatorTraits to provide the field type calculators. Builders and
 * evaluators acting on several index ranges at once (e.g. on a multipatch geometry)
 * provide their own field type calculators (e.g. MultipatchField).
 *
 * @tparam MemorySpace The memory space where the data is saved.
 */
template <class MemorySpace>
struct InterpolationFieldTypes
{
    /// @brief A modifiable field of ElementType defined on IdxRangeType.
    template <class ElementType, class IdxRangeType>
    using field_type = Field<ElementType, IdxRangeType, MemorySpace>;

    /// @brief A constant field of ElementType defined on IdxRangeType.
    template <class ElementType, class IdxRangeType>
    using const_field_type = ConstField<ElementType, IdxRangeType, MemorySpace>;

    /// @brief A modifiable vector field of ElementType defined on IdxRangeType.
    template <class ElementType, class IdxRangeType, class VectorIndexSetType>
    using vector_field_type
            = VectorField<ElementType, IdxRangeType, VectorIndexSetType, MemorySpace>;

    /// @brief A constant vector field of ElementType defined on IdxRangeType.
    template <class ElementType, class IdxRangeType, class VectorIndexSetType>
    using const_vector_field_type
            = VectorConstField<ElementType, IdxRangeType, VectorIndexSetType, MemorySpace>;
};

} // namespace detail
