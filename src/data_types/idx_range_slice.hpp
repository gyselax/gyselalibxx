// SPDX-License-Identifier: MIT
#pragma once
#include <map>
#include <vector>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


template <class...>
class IdxRangeSlice;

template <class Grid1D>
struct IdxRangeSliceIterator;

template <class T>
struct is_subidx_range_collection : std::false_type
{
};

template <class... Tags>
struct is_subidx_range_collection<IdxRangeSlice<Tags...>> : std::true_type
{
};

/// @brief A templated value to determine if a type is a IdxRangeSlice.
template <class T>
inline constexpr bool is_subidx_range_collection_v = is_subidx_range_collection<T>::value;

/**
 * @brief A class which describes a collection of equally spaced Idxs which form a index range.
 *
 * This class should eventually be replaced by a DDC functionality when this becomes available.
 */
template <class... Dims>
class IdxRangeSlice
{
    static_assert(sizeof...(Dims) > 0);

private:
    /// @brief A Type Sequence indicating the order of the dimensions
    using tags_seq = ddc::detail::TypeSeq<Dims...>;

    Idx<Dims...> m_front;
    IdxStep<Dims...> m_size;
    IdxStep<Dims...> m_stride;

    template <class...>
    friend class IdxRangeSlice;

public:
    /// @brief Default constructor for IdxRangeSlice creating an empty index range.
    IdxRangeSlice() = default;

    /**
     * @brief Build a IdxRangeSlice from vectors of valid Idxs in each dimension.
     *
     * @param front The Idx describing the first point in the index range.
     * @param size The number of elements along each direction.
     * @param stride The IdxStep distance between subsequent elements of the index range.
     */
    KOKKOS_FUNCTION IdxRangeSlice(
            Idx<Dims...> front,
            IdxStep<Dims...> size,
            IdxStep<Dims...> stride)
        : m_front(front)
        , m_size(size)
        , m_stride(stride)
    {
    }

    /**
     * @brief Build a IdxRangeSlice from a set of 1D IdxRangeSlices.
     *
     * @param valid_indices The IdxRangeSlices which comprise this IdxRangeSlice.
     */
    template <
            class... DDoms,
            class = std::enable_if_t<(is_subidx_range_collection_v<DDoms> && ...)>>
    KOKKOS_FUNCTION IdxRangeSlice(DDoms const&... valid_indices)
        : m_front(valid_indices.front()...)
        , m_size(valid_indices.extents()...)
        , m_stride(valid_indices.strides()...)
    {
    }

    /**
     * @brief Check if the specified Idx is found in this index range slice.
     *
     * @param elem The element which may or may not be in this index range slice.
     *
     * @returns bool True if the element is found in this index range slice, False otherwise.
     */
    template <class... DDims>
    KOKKOS_FUNCTION bool contains(Idx<DDims...> elem) const
    {
        static_assert(
                (ddc::in_tags_v<DDims, tags_seq> && ...),
                "Requested Tag absent from IdxRangeSlice");

        IdxStep<DDims...> distance(elem - ddc::select<DDims...>(m_front));
        return (((ddc::get<DDims>(distance) % ddc::get<DDims>(m_stride)) == 0) && ...);
    }

    /**
     * @brief Check if all elements of the specified IdxRange are found in this index range slice.
     *
     * @param idx_range The index range which may or may not be in this index range slice.
     *
     * @returns bool True if the all elements of the specified IdxRange are found in this index range slice.
     *               False otherwise.
     */
    template <class... DDims>
    KOKKOS_FUNCTION bool contains(IdxRange<DDims...> idx_range) const
    {
        static_assert(
                (ddc::in_tags_v<DDims, tags_seq> && ...),
                "Requested Tag absent from IdxRangeSlice");
        return contains(idx_range.front())
               && (((idx_range.template extent<DDims>().value() == 1)
                    || (ddc::select<DDims>(m_stride) == 1))
                   && ...);
    }

    /**
     * @brief Get the index of the Idx within the index range slice.
     * This function is particularly useful to index an mdspan over the index range slice.
     *
     * @param elem A 1D Idx which is inside the index range slice.
     *
     * @returns The index of the element.
     */
    template <class Dim, class = std::enable_if_t<ddc::in_tags_v<Dim, tags_seq>>>
    KOKKOS_FUNCTION std::size_t get_index(Idx<Dim> elem) const
    {
        assert(contains(elem));
        return (elem - ddc::select<Dim>(m_front)).value() / ddc::get<Dim>(m_stride);
    }

    /**
     * @brief Get the size of the index range slice in each dimension.
     *
     * @returns A IdxStep describing the size of the index range slice in each dimension.
     */
    KOKKOS_FUNCTION constexpr IdxStep<Dims...> extents() const noexcept
    {
        return m_size;
    }

    /**
     * @brief Get the size of the index range slice in the specified dimension.
     *
     * @tparam QueryDDim The dimension in which the size is requested.
     *
     * @returns A IdxStep describing the size of the index range slice in the specified dimension.
     */
    template <class QueryDDim>
    KOKKOS_FUNCTION constexpr IdxStep<QueryDDim> extent() const noexcept
    {
        return ddc::select<QueryDDim>(m_size);
    }

    /**
     * @brief Get the total number of elements in the index range slice.
     *
     * @return The total number of elements in the index range slice.
     */
    KOKKOS_FUNCTION constexpr std::size_t size() const
    {
        return (1ul * ... * (ddc::get<Dims>(m_size)));
    }

    /**
     * @brief Get the first element in the index range slice.
     *
     * @return The first element in the index range slice.
     */
    KOKKOS_FUNCTION constexpr Idx<Dims...> front() const noexcept
    {
        return m_front;
    }

    /**
     * @brief Get the last element in the index range slice.
     *
     * @return The last element in the index range slice.
     */
    KOKKOS_FUNCTION constexpr Idx<Dims...> back() const noexcept
    {
        return m_front + m_stride * (m_size - 1);
    }

    /**
     * @brief Get the stride from one element of the index range slice to another.
     *
     * @tparam QueryDim The dimension being queried
     *
     * @return A 1D IdxStep describing a stride from one element to another.
     */
    template <class QueryDim>
    KOKKOS_FUNCTION constexpr IdxStep<QueryDim> stride() const noexcept
    {
        return ddc::select<QueryDim>(m_stride);
    }

    /**
     * @brief Get the strides from one element of the index range slice to another.
     *
     * @return A IdxStep describing the strides from one element to another.
     */
    KOKKOS_FUNCTION constexpr IdxStep<Dims...> strides() const noexcept
    {
        return m_stride;
    }

    /**
     * @brief Get the iterator to the first element of the IdxRangeSlice.
     * The elements are pairs of Idxs and indices.
     *
     * @returns The iterator to the first element of the IdxRangeSlice.
     */
    KOKKOS_FUNCTION auto begin() const
    {
        static_assert(sizeof...(Dims) == 1);
        return IdxRangeSliceIterator(m_front, m_stride);
    }

    /**
     * @brief Get the iterator to the end of the IdxRangeSlice.
     * The elements are pairs of Idxs and indices.
     *
     * @returns The iterator to the end of the IdxRangeSlice.
     */
    KOKKOS_FUNCTION auto end() const
    {
        static_assert(sizeof...(Dims) == 1);
        return IdxRangeSliceIterator(m_front + m_stride * m_size, m_stride);
    }
};

/**
 * @brief An iterator type for the IdxRangeSlice.
 */
template <class Grid1D>
struct IdxRangeSliceIterator
{
private:
    Idx<Grid1D> m_value;
    IdxStep<Grid1D> m_stride;

public:
    /// The type of iterator
    using iterator_category = std::random_access_iterator_tag;

    /// The type of the values stored in the iterator
    using value_type = Idx<Grid1D>;

    /// The type of the stride between values
    using stride_type = IdxStep<Grid1D>;

    /// The type that can be used to increment the iterator
    using difference_type = std::ptrdiff_t;

    KOKKOS_DEFAULTED_FUNCTION IdxRangeSliceIterator() = default;

    /**
     * @brief Build an iterator using the current value and the distance to the following element
     *
     * @param value The value of the discrete sub-index range element.
     * @param stride The stride between consectuive sub-index range elements.
     */
    KOKKOS_FUNCTION constexpr explicit IdxRangeSliceIterator(
            Idx<Grid1D> value,
            IdxStep<Grid1D> stride)
        : m_value(value)
        , m_stride(stride)
    {
    }

    /**
     * @brief Get the value referred to by the iterator.
     * @return The value referred to by the iterator.
     */
    KOKKOS_FUNCTION constexpr Idx<Grid1D> operator*() const noexcept
    {
        return m_value;
    }

    /**
     * @brief The prefix increment operator.
     * @return A reference to the current incremented iterator.
     */
    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator& operator++()
    {
        m_value = m_value + m_stride;
        return *this;
    }

    /**
     * @brief The postfix increment operator.
     * @return A reference to the non-incremented iterator.
     */
    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator operator++(int)
    {
        auto tmp = *this;
        ++*this;
        return tmp;
    }

    /**
     * @brief The prefix decrement operator.
     * @return A reference to the current decremented iterator.
     */
    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator& operator--()
    {
        m_value = m_value - m_stride;
        return *this;
    }

    /**
     * @brief The postfix decrement operator.
     * @return A reference to the non-decremented iterator.
     */
    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator operator--(int)
    {
        auto tmp = *this;
        --*this;
        return tmp;
    }

    /**
     * @brief Increment the current iterator by n elements.
     * @param n The number of strides by which the iterator should be incremented.
     * @return A reference to the current incremented iterator.
     */
    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator& operator+=(difference_type n)
    {
        m_value = m_value + n * m_stride;
        return *this;
    }

    /**
     * @brief Decrement the current iterator by n elements.
     * @param n The number of strides by which the iterator should be decremented.
     * @return A reference to the current decremented iterator.
     */
    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator& operator-=(difference_type n)
    {
        m_value = m_value - n * m_stride;
        return *this;
    }

    /**
     * @brief Access the n-th following element.
     * @param n The index of the element.
     * @return The Idx n-th position in the sub-index range following the value
     *          indicated by this iterator.
     */
    KOKKOS_FUNCTION constexpr Idx<Grid1D> operator[](difference_type n) const
    {
        return m_value + n * m_stride;
    }

    /// Compare iterator equality
    friend KOKKOS_FUNCTION constexpr bool operator==(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return xx.m_value == yy.m_value;
    }

    /// Compare iterator non-equality
    friend KOKKOS_FUNCTION constexpr bool operator!=(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return xx.m_value != yy.m_value;
    }

    /// Compare the order of iterators
    friend KOKKOS_FUNCTION constexpr bool operator<(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return xx.m_value < yy.m_value;
    }

    /// Compare the order of iterators
    friend KOKKOS_FUNCTION constexpr bool operator>(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return yy < xx;
    }

    /// Compare the order of iterators
    friend KOKKOS_FUNCTION constexpr bool operator<=(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return !(yy < xx);
    }

    /// Compare the order of iterators
    friend KOKKOS_FUNCTION constexpr bool operator>=(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return !(xx < yy);
    }

    /// Increment an iterator by n elements.
    friend KOKKOS_FUNCTION constexpr IdxRangeSliceIterator operator+(
            IdxRangeSliceIterator i,
            difference_type n)
    {
        return i += n;
    }

    /// Increment an iterator by n elements.
    friend KOKKOS_FUNCTION constexpr IdxRangeSliceIterator operator+(
            difference_type n,
            IdxRangeSliceIterator i)
    {
        return i += n;
    }

    /// Decrement an iterator by n elements.
    friend KOKKOS_FUNCTION constexpr IdxRangeSliceIterator operator-(
            IdxRangeSliceIterator i,
            difference_type n)
    {
        return i -= n;
    }

    /// Decrement an iterator by n elements.
    friend KOKKOS_FUNCTION constexpr difference_type operator-(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        assert(xx.m_stride == yy.m_stride);
        return (yy.m_value > xx.m_value)
                       ? ((-static_cast<difference_type>(yy.m_value - xx.m_value))
                          / ddc::get<Grid1D>(xx.m_stride))
                       : ((xx.m_value - yy.m_value) / ddc::get<Grid1D>(xx.m_stride));
    }
};

/// A class to create a IdxRangeSlice type from a TypeSeq.
template <class T>
struct IdxRangeToSlice;

template <class... Dims>
struct IdxRangeToSlice<ddc::detail::TypeSeq<Dims...>>
{
    using value = IdxRangeSlice<Dims...>;
};

/// @brief A templated type for creating a IdxRangeSlice from a TypeSeq.
template <class T>
using to_subidx_range_collection = typename IdxRangeToSlice<T>::value;
