// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "tensor.hpp"

template <class FieldType, class NDTypeSeq>
class VectorFieldCommon;

template <class T>
inline constexpr bool enable_vector_field = false;

template <class T>
inline constexpr bool enable_borrowed_vector_field = false;

template <class T>
inline constexpr bool is_vector_field_v
        = enable_vector_field<std::remove_const_t<std::remove_reference_t<T>>>;

template <class T>
inline constexpr bool is_borrowed_vector_field_v
        = is_vector_field_v<
                  T> && (std::is_lvalue_reference_v<T> || enable_borrowed_vector_field<std::remove_cv_t<std::remove_reference_t<T>>>);

namespace ddcHelper {

template <
        class FieldDst,
        class FieldSrc,
        std::enable_if_t<
                is_borrowed_vector_field_v<FieldDst> && is_borrowed_vector_field_v<FieldSrc>,
                bool> = true>
auto deepcopy(FieldDst&& dst, FieldSrc&& src)
{
    static_assert(std::is_same_v<
                  typename std::remove_reference_t<FieldDst>::NDTypeTag,
                  typename std::remove_reference_t<FieldSrc>::NDTypeTag>);

    assert(dst.idx_range().extents() == src.idx_range().extents());

    dst.deepcopy(src);
    return get_field(dst);
}

template <
        class QueryTag,
        class VectorFieldType,
        std::enable_if_t<is_vector_field_v<VectorFieldType>, bool> = true>
inline constexpr auto get(VectorFieldType const& field) noexcept
{
    static_assert(is_vector_field_v<VectorFieldType>);
    return field.template get<QueryTag>();
}

template <
        class QueryTag,
        class VectorFieldType,
        std::enable_if_t<is_vector_field_v<VectorFieldType>, bool> = true>
inline constexpr auto get(VectorFieldType& field) noexcept
{
    static_assert(is_vector_field_v<VectorFieldType>);
    return field.template get<QueryTag>();
}


} // namespace ddcHelper

/**
 * @brief A base class containing the functionalities common to VectorField and VectorFieldMem.
 *
 * The vector field is stored as one scalar field per component.
 *
 * @tparam FieldType The type of the scalar field used to store each component.
 * @tparam DDims The tags describing the dimensions which make up the elements of the vector.
 */
template <class FieldType, class... DDims>
class VectorFieldCommon<FieldType, ddc::TypeSeq<DDims...>>
{
    static_assert(ddc::is_chunk_v<FieldType>);
    using data_type = typename FieldType::element_type;

public:
    /// @brief The type of the elements in the fields.
    using element_type = Vector<std::remove_const_t<data_type>, DDims...>;

    /**
     * @brief The IdxRange on which the fields in this object are defined.
     *
     * This is a DDC keyword used to make this class interchangeable with Field.
     * In DDC IdxRange types are referred to as DiscreteDomain types.
     */
    using discrete_domain_type = typename FieldType::discrete_domain_type;
    /// @brief The IdxRange on which the fields in this object are defined.
    using index_range_type = discrete_domain_type;

    /// @brief The type of the element returned when indexing this field (e.g. a tuple of doubles).
    using element_ref_type = Vector<data_type&, DDims...>;

    /// @brief The tags describing the dimensions which make up the elements of the vector.
    using NDTypeTag = ddc::TypeSeq<DDims...>;

    /// @brief The type of a modifiable span of one of these field. This is a DDC keyword used to make this class interchangeable with Field.
    using chunk_span_type = typename FieldType::span_type;
    /// @brief The type of a constant view of one of these field. This is a DDC keyword used to make this class interchangeable with Field.
    using chunk_view_type = typename FieldType::view_type;

protected:
    /// @brief The number of components of the vector.
    static constexpr std::size_t NDims = sizeof...(DDims);

    /// @brief The scalar fields storing each component of the vector field.
    std::array<FieldType, NDims> m_values;

protected:
    /// Empty Field
    KOKKOS_DEFAULTED_FUNCTION VectorFieldCommon() = default;

    /// Construct a Field on an index range with uninitialized values
    template <
            class... FieldTypes,
            class = std::enable_if_t<std::conjunction_v<std::is_same<FieldTypes, FieldType>...>>,
            std::enable_if_t<
                    !std::conjunction_v<
                            std::bool_constant<ddc::is_borrowed_chunk_v<FieldTypes>>...>,
                    int> = 0>
    explicit VectorFieldCommon(FieldTypes&&... fields) : m_values {std::move(fields)...}
    {
    }

    /// Construct a Field on an index range with uninitialized values
    template <
            class... FieldTypes,
            class = std::enable_if_t<std::conjunction_v<std::is_same<FieldTypes, FieldType>...>>,
            std::enable_if_t<
                    std::conjunction_v<std::bool_constant<ddc::is_borrowed_chunk_v<FieldTypes>>...>,
                    int> = 0>
    KOKKOS_FUNCTION explicit VectorFieldCommon(FieldTypes&&... fields)
        : m_values {std::forward<FieldTypes>(fields)...}
    {
    }

public:
    /**
     * @brief Provide access to the index range on which this field is defined.
     * @return The index range on which this field is defined.
     */
    constexpr index_range_type idx_range() const noexcept
    {
        return m_values[0].domain();
    }

    /**
     * @brief Provide access to the sub-set of the index range on which this field is defined.
     * @tparam QueryDDims The dimensions of interest.
     * @return The index range on which this field is defined.
     */
    template <class... QueryDDims>
    constexpr IdxRange<QueryDDims...> idx_range() const noexcept
    {
        return ddc::select<QueryDDims...>(idx_range());
    }

    /**
     * @brief Get the number of dimensions of the index range on which the field is defined.
     * @return The rank of the index range.
     */
    static constexpr int rank() noexcept
    {
        return FieldType::rank();
    }

    /**
     * @brief Copy the values of another vector field into this vector field.
     * @param[in] src The vector field whose values should be copied.
     */
    template <class FieldSrc>
    void deepcopy(FieldSrc const& src)
    {
        ((ddc::parallel_deepcopy(
                 get_field(m_values[ddc::type_seq_rank_v<DDims, NDTypeTag>]),
                 ddcHelper::get<DDims>(src))),
         ...);
    }
};
