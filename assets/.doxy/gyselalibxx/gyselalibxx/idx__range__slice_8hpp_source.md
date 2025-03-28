

# File idx\_range\_slice.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**idx\_range\_slice.hpp**](idx__range__slice_8hpp.md)

[Go to the documentation of this file](idx__range__slice_8hpp.md)


```C++
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

template <class T>
inline constexpr bool is_subidx_range_collection_v = is_subidx_range_collection<T>::value;

template <class... Dims>
class IdxRangeSlice
{
    static_assert(sizeof...(Dims) > 0);

private:
    using tags_seq = ddc::detail::TypeSeq<Dims...>;

    Idx<Dims...> m_front;
    IdxStep<Dims...> m_size;
    IdxStep<Dims...> m_stride;

    template <class...>
    friend class IdxRangeSlice;

public:
    IdxRangeSlice() = default;

    KOKKOS_FUNCTION IdxRangeSlice(
            Idx<Dims...> front,
            IdxStep<Dims...> size,
            IdxStep<Dims...> stride)
        : m_front(front)
        , m_size(size)
        , m_stride(stride)
    {
    }

    template <
            class... DDoms,
            class = std::enable_if_t<(is_subidx_range_collection_v<DDoms> && ...)>>
    KOKKOS_FUNCTION IdxRangeSlice(DDoms const&... valid_indices)
        : m_front(valid_indices.front()...)
        , m_size(valid_indices.extents()...)
        , m_stride(valid_indices.strides()...)
    {
    }

    template <class... DDims>
    KOKKOS_FUNCTION bool contains(Idx<DDims...> elem) const
    {
        static_assert(
                (ddc::in_tags_v<DDims, tags_seq> && ...),
                "Requested Tag absent from IdxRangeSlice");

        IdxStep<DDims...> distance(elem - ddc::select<DDims...>(m_front));
        return (((ddc::get<DDims>(distance) % ddc::get<DDims>(m_stride)) == 0) && ...);
    }

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

    template <class Dim, class = std::enable_if_t<ddc::in_tags_v<Dim, tags_seq>>>
    KOKKOS_FUNCTION std::size_t get_index(Idx<Dim> elem) const
    {
        assert(contains(elem));
        return (elem - ddc::select<Dim>(m_front)).value() / ddc::get<Dim>(m_stride);
    }

    KOKKOS_FUNCTION constexpr IdxStep<Dims...> extents() const noexcept
    {
        return m_size;
    }

    template <class QueryDDim>
    KOKKOS_FUNCTION constexpr IdxStep<QueryDDim> extent() const noexcept
    {
        return ddc::select<QueryDDim>(m_size);
    }

    KOKKOS_FUNCTION constexpr std::size_t size() const
    {
        return (1ul * ... * (ddc::get<Dims>(m_size)));
    }

    KOKKOS_FUNCTION constexpr Idx<Dims...> front() const noexcept
    {
        return m_front;
    }

    KOKKOS_FUNCTION constexpr Idx<Dims...> back() const noexcept
    {
        return m_front + m_stride * (m_size - 1);
    }

    template <class QueryDim>
    KOKKOS_FUNCTION constexpr IdxStep<QueryDim> stride() const noexcept
    {
        return ddc::select<QueryDim>(m_stride);
    }

    KOKKOS_FUNCTION constexpr IdxStep<Dims...> strides() const noexcept
    {
        return m_stride;
    }

    KOKKOS_FUNCTION auto begin() const
    {
        static_assert(sizeof...(Dims) == 1);
        return IdxRangeSliceIterator(m_front, m_stride);
    }

    KOKKOS_FUNCTION auto end() const
    {
        static_assert(sizeof...(Dims) == 1);
        return IdxRangeSliceIterator(m_front + m_stride * m_size, m_stride);
    }
};

template <class Grid1D>
struct IdxRangeSliceIterator
{
private:
    Idx<Grid1D> m_value;
    IdxStep<Grid1D> m_stride;

public:
    using iterator_category = std::random_access_iterator_tag;

    using value_type = Idx<Grid1D>;

    using stride_type = IdxStep<Grid1D>;

    using difference_type = std::ptrdiff_t;

    KOKKOS_DEFAULTED_FUNCTION IdxRangeSliceIterator() = default;

    KOKKOS_FUNCTION constexpr explicit IdxRangeSliceIterator(
            Idx<Grid1D> value,
            IdxStep<Grid1D> stride)
        : m_value(value)
        , m_stride(stride)
    {
    }

    KOKKOS_FUNCTION constexpr Idx<Grid1D> operator*() const noexcept
    {
        return m_value;
    }

    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator& operator++()
    {
        m_value = m_value + m_stride;
        return *this;
    }

    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator operator++(int)
    {
        auto tmp = *this;
        ++*this;
        return tmp;
    }

    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator& operator--()
    {
        m_value = m_value - m_stride;
        return *this;
    }

    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator operator--(int)
    {
        auto tmp = *this;
        --*this;
        return tmp;
    }

    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator& operator+=(difference_type n)
    {
        m_value = m_value + n * m_stride;
        return *this;
    }

    KOKKOS_FUNCTION constexpr IdxRangeSliceIterator& operator-=(difference_type n)
    {
        m_value = m_value - n * m_stride;
        return *this;
    }

    KOKKOS_FUNCTION constexpr Idx<Grid1D> operator[](difference_type n) const
    {
        return m_value + n * m_stride;
    }

    friend KOKKOS_FUNCTION constexpr bool operator==(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return xx.m_value == yy.m_value;
    }

    friend KOKKOS_FUNCTION constexpr bool operator!=(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return xx.m_value != yy.m_value;
    }

    friend KOKKOS_FUNCTION constexpr bool operator<(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return xx.m_value < yy.m_value;
    }

    friend KOKKOS_FUNCTION constexpr bool operator>(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return yy < xx;
    }

    friend KOKKOS_FUNCTION constexpr bool operator<=(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return !(yy < xx);
    }

    friend KOKKOS_FUNCTION constexpr bool operator>=(
            IdxRangeSliceIterator const& xx,
            IdxRangeSliceIterator const& yy)
    {
        return !(xx < yy);
    }

    friend KOKKOS_FUNCTION constexpr IdxRangeSliceIterator operator+(
            IdxRangeSliceIterator i,
            difference_type n)
    {
        return i += n;
    }

    friend KOKKOS_FUNCTION constexpr IdxRangeSliceIterator operator+(
            difference_type n,
            IdxRangeSliceIterator i)
    {
        return i += n;
    }

    friend KOKKOS_FUNCTION constexpr IdxRangeSliceIterator operator-(
            IdxRangeSliceIterator i,
            difference_type n)
    {
        return i -= n;
    }

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

template <class T>
struct IdxRangeToSlice;

template <class... Dims>
struct IdxRangeToSlice<ddc::detail::TypeSeq<Dims...>>
{
    using value = IdxRangeSlice<Dims...>;
};

template <class T>
using to_subidx_range_collection = typename IdxRangeToSlice<T>::value;
```


