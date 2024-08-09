// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

template <class, class>
class VectorFieldCommon;

template <class T>
inline constexpr bool enable_field = false;

template <class T>
inline constexpr bool enable_borrowed_field = false;

template <class T>
inline constexpr bool is_field_v = enable_field<std::remove_const_t<std::remove_reference_t<T>>>;

template <class T>
inline constexpr bool is_borrowed_field_v
        = is_field_v<
                  T> && (std::is_lvalue_reference_v<T> || enable_borrowed_field<std::remove_cv_t<std::remove_reference_t<T>>>);

namespace ddcHelper {

template <class... QueryDDims, class FieldType>
auto get_domain(FieldType const& field) noexcept
{
    static_assert(is_field_v<FieldType>, "Not a field span type");
    return field.template domain<QueryDDims...>();
}

template <
        class FieldDst,
        class FieldSrc,
        std::enable_if_t<
                is_borrowed_field_v<FieldDst> && is_borrowed_field_v<FieldSrc>,
                bool> = true>
auto deepcopy(FieldDst&& dst, FieldSrc&& src)
{
    static_assert(std::is_same_v<
                  typename std::remove_reference_t<FieldDst>::NDTypeTag,
                  typename std::remove_reference_t<FieldSrc>::NDTypeTag>);

    assert(dst.domain().extents() == src.domain().extents());

    dst.deepcopy(src);
    return dst.span_view();
}

template <class QueryTag, class ChunkType, class NDTypeTag>
inline constexpr typename ChunkType::span_type get(
        VectorFieldCommon<ChunkType, NDTypeTag>& field) noexcept
{
    return field.template get<QueryTag>();
}

template <class QueryTag, class ChunkType, class NDTypeTag>
inline constexpr typename ChunkType::view_type get(
        VectorFieldCommon<ChunkType, NDTypeTag> const& field) noexcept
{
    return field.template get<QueryTag>();
}


} // namespace ddcHelper

template <class ChunkType, class... DDims>
class VectorFieldCommon<ChunkType, ddc::detail::TypeSeq<DDims...>>
{
    static_assert(ddc::is_chunk_v<ChunkType>);
    using data_type = typename ChunkType::element_type;

public:
    using discrete_domain_type = typename ChunkType::discrete_domain_type;
    using element_type = typename ddc::detail::TaggedVector<data_type, DDims...>;
    using element_ref_type = typename ddc::detail::TaggedVector<data_type&, DDims...>;
    using NDTypeTag = ddc::detail::TypeSeq<DDims...>;

    using chunk_span_type = typename ChunkType::span_type;
    using chunk_view_type = typename ChunkType::view_type;

protected:
    static constexpr std::size_t NDims = sizeof...(DDims);

    std::array<ChunkType, NDims> m_values;

protected:
    /// Empty Chunk
    KOKKOS_DEFAULTED_FUNCTION VectorFieldCommon() = default;

    /// Construct a Chunk on a domain with uninitialized values
    template <
            class... Chunks,
            class = std::enable_if_t<std::conjunction_v<std::is_same<Chunks, ChunkType>...>>,
            std::enable_if_t<
                    !std::conjunction_v<std::bool_constant<ddc::is_borrowed_chunk_v<Chunks>>...>,
                    int> = 0>
    explicit VectorFieldCommon(Chunks&&... chunks) : m_values {std::move(chunks)...}
    {
    }

    /// Construct a Chunk on a domain with uninitialized values
    template <
            class... Chunks,
            class = std::enable_if_t<std::conjunction_v<std::is_same<Chunks, ChunkType>...>>,
            std::enable_if_t<
                    std::conjunction_v<std::bool_constant<ddc::is_borrowed_chunk_v<Chunks>>...>,
                    int> = 0>
    KOKKOS_FUNCTION explicit VectorFieldCommon(Chunks&&... chunks)
        : m_values {std::forward<Chunks>(chunks)...}
    {
    }

public:
    constexpr discrete_domain_type domain() const noexcept
    {
        return m_values[0].domain();
    }

    /** Provide access to the domain on which this chunk is defined
     * @return the domain on which this chunk is defined
     */
    template <class... QueryDDims>
    constexpr ddc::DiscreteDomain<QueryDDims...> domain() const noexcept
    {
        return ddc::select<QueryDDims...>(domain());
    }

    static constexpr int rank() noexcept
    {
        return ChunkType::rank();
    }

    template <class FieldSrc>
    void deepcopy(FieldSrc const& src)
    {
        ((ddc::parallel_deepcopy(this->get<DDims>(), src.template get<DDims>())), ...);
    }

    template <class QueryTag>
    inline constexpr chunk_span_type get() noexcept
    {
        static_assert(
                ddc::in_tags_v<QueryTag, NDTypeTag>,
                "requested Tag absent from TaggedVector");
        return m_values[ddc::type_seq_rank_v<QueryTag, NDTypeTag>].span_view();
    }

    template <class QueryTag>
    inline constexpr chunk_view_type get() const noexcept
    {
        static_assert(
                ddc::in_tags_v<QueryTag, NDTypeTag>,
                "requested Tag absent from TaggedVector");
        return m_values[ddc::type_seq_rank_v<QueryTag, NDTypeTag>].span_cview();
    }
};
