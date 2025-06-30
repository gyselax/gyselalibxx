// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "tensor.hpp"
#include "vector_field_common.hpp"

template <
        class ElementType,
        class IdxRangeType,
        class VectorIndexSetType,
        class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
class VectorFieldMem
    : private FieldMem<
              ElementType,
              ddc::cartesian_prod_t<
                      IdxRange<detail::VecIdxSetWrapper<VectorIndexSetType>>,
                      IdxRangeType>,
              MemSpace>
{
    static_assert(is_vector_index_set_v<VectorIndexSetType>);

private:
    using base_type = FieldMem<
            ElementType,
            ddc::cartesian_prod_t<
                    IdxRange<detail::VecIdxSetWrapper<VectorIndexSetType>>,
                    IdxRangeType>,
            MemSpace>;

public:
    using idx_range_type = IdxRangeType;
    using Allocator = ddc::KokkosAllocator<ElementType, MemSpace>;
    /// @brief The type of the elements in the fields.
    using element_type = Tensor<std::remove_const_t<ElementType>, VectorIndexSetType>;
    using base_type::layout_type;

private:
    using internal_idx_range_type = typename base_type::discrete_domain_type;

public:
    /// Construct a labeled VectorFieldMem on a domain with uninitialized values
    explicit VectorFieldMem(
            std::string const& label,
            idx_range_type const& domain,
            Allocator allocator = Allocator())
        : base_type(
                label,
                internal_idx_range_type(detail::make_idx_range<VectorIndexSetType>(), domain),
                allocator)
    {
    }

    /// Construct a VectorFieldMem on a domain with uninitialized values
    explicit VectorFieldMem(idx_range_type const& domain, Allocator allocator = Allocator())
        : VectorFieldMem("no-label", domain, std::move(allocator))
    {
    }

    idx_range_type idx_range() const
    {
        return idx_range_type(this->domain());
    }

    /** Element access using a multi-dimensional Idx
     * @param delems discrete coordinates
     * @return copy of this element
     */
    template <class... ODDims, typename T, T... ints>
    element_type const operator()(Idx<ODDims...> const& delems, std::integer_sequence<T, ints...>)
            const noexcept
    {
        return element_type((base_type::m_values[ints](delems))...);
    }



    using base_type::span_view;
    using base_type::span_cview;
};


template <class ElementType, class IdxRangeType, class DimSeq, class MemSpace>
inline constexpr bool
        enable_vector_field<VectorFieldMem<ElementType, IdxRangeType, DimSeq, MemSpace>> = true;

template <class ElementType, class IdxRangeType, class DimSeq, class Allocator>
inline constexpr bool enable_data_access_methods<
        VectorFieldMem<ElementType, IdxRangeType, DimSeq, Allocator>> = true;

template <class ElementType, class IdxRangeType, class DimSeq, class Allocator>
inline constexpr bool
        enable_mem_type<VectorFieldMem<ElementType, IdxRangeType, DimSeq, Allocator>> = true;

namespace detail {
/**
 * @brief Set a `VectorFieldMem` on a given NewMemorySpace.
 * @tparam NewMemorySpace The new memory space.
 * @tparam ElementType Type of the elements in the ddc::Chunk of the VectorFieldMem.
 * @tparam SupportType Type of the domain of the ddc::Chunk in the VectorFieldMem.
 * @tparam VectorIndexSetType VectorIndexSet object storing the dimensions along which the VectorFieldMem is defined.
 *               The dimensions refer to the dimensions of the arrival domain of the VectorFieldMem.
 * @tparam MemSpace The old memory space.
 * @see VectorFieldMem
 */
template <
        class NewMemorySpace,
        class ElementType,
        class SupportType,
        class VectorIndexSetType,
        class MemSpace>
struct OnMemorySpace<
        NewMemorySpace,
        VectorFieldMem<ElementType, SupportType, VectorIndexSetType, MemSpace>>
{
    using type = VectorFieldMem<ElementType, SupportType, VectorIndexSetType, NewMemorySpace>;
};
} // namespace detail

namespace ddcHelper {

template <
        class QueryTag,
        class ElementType,
        class SupportType,
        class VectorIndexSetType,
        class MemorySpace>
inline constexpr auto get(
        VectorFieldMem<ElementType, SupportType, VectorIndexSetType, MemorySpace> const&
                field) noexcept
{
    Idx<detail::VecIdxSetWrapper<VectorIndexSetType>> idx(
            ddc::type_seq_rank_v<QueryTag, VectorIndexSetType>);
    return field[idx];
}

template <
        class QueryTag,
        class ElementType,
        class SupportType,
        class VectorIndexSetType,
        class MemorySpace>
inline constexpr auto get(
        VectorFieldMem<ElementType, SupportType, VectorIndexSetType, MemorySpace>& field) noexcept
{
    Idx<detail::VecIdxSetWrapper<VectorIndexSetType>> idx(
            ddc::type_seq_rank_v<QueryTag, VectorIndexSetType>);
    return field[idx];
}

} // namespace ddcHelper
