// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "vector_field_common.hpp"

template <
        class ElementType,
        class IdxRangeType,
        class VectorIndexSetType,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using VectorField
        = Field<ElementType,
                ddc::cartesian_prod_t<IdxRange<detail::VecIdxSetWrapper<VectorIndexSetType>>, IdxRangeType>,
                MemorySpace,
                LayoutStridedPolicy>;

namespace ddcHelper {

template <
        class QueryTag,
        class ElementType,
        class VectorIndexSetType,
        class... Grid,
        class MemorySpace,
        class LayoutStridedPolicy>
inline constexpr auto get(
        Field<ElementType,
              IdxRange<detail::VecIdxSetWrapper<VectorIndexSetType>, Grid...>,
              MemorySpace,
              LayoutStridedPolicy> const& field) noexcept
{
    Idx<detail::VecIdxSetWrapper<VectorIndexSetType>> idx(
            ddc::type_seq_rank_v<QueryTag, VectorIndexSetType>);
    return field[idx];
}

} // namespace ddcHelper
