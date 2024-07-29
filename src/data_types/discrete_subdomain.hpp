#pragma once
#include <map>
#include <vector>

#include <ddc/ddc.hpp>


template <class...>
class DiscreteSubDomain;

template <class Grid1D>
struct DiscreteSubDomainIterator;

template <class T>
struct is_subdomain_collection : std::false_type
{
};

template <class... Tags>
struct is_subdomain_collection<DiscreteSubDomain<Tags...>> : std::true_type
{
};

/// @brief A templated value to determine if a type is a DiscreteSubDomain.
template <class T>
inline constexpr bool is_subdomain_collection_v = is_subdomain_collection<T>::value;

/**
 * @brief A class which describes a collection of ddc::DiscreteElements which form a subdomain.
 *
 * This class should eventually be replaced by a DDC functionality when this becomes available.
 */
template <class... Dims>
class DiscreteSubDomain
{
    static_assert(sizeof...(Dims) > 0);

private:
    /// @brief A Type Sequence indicating the order of the dimensions
    using tags_seq = ddc::detail::TypeSeq<Dims...>;

    ddc::DiscreteElement<Dims...> m_front;
    ddc::DiscreteVector<Dims...> m_size;
    ddc::DiscreteVector<Dims...> m_stride;

    template <class...>
    friend class DiscreteSubDomain;

public:
    /// @brief Default constructor for DiscreteSubDomain creating an empty domain.
    DiscreteSubDomain() = default;

    /**
     * @brief Build a DiscreteSubDomain from vectors of valid DiscreteElements in each dimension.
     *
     * @param front The DiscreteElement describing the first point in the domain.
     * @param size The number of elements along each direction.
     * @param stride The DiscreteVector distance between subsequent elements of the domain.
     */
    KOKKOS_FUNCTION DiscreteSubDomain(
            ddc::DiscreteElement<Dims...> front,
            ddc::DiscreteVector<Dims...> size,
            ddc::DiscreteVector<Dims...> stride)
        : m_front(front)
        , m_size(size)
        , m_stride(stride)
    {
    }

    /**
     * @brief Build a DiscreteSubDomain from a set of 1D DiscreteSubDomains.
     *
     * @param valid_indices The DiscreteSubDomains which comprise this DiscreteSubDomain.
     */
    template <class... DDoms, class = std::enable_if_t<(is_subdomain_collection_v<DDoms> && ...)>>
    KOKKOS_FUNCTION DiscreteSubDomain(DDoms const&... valid_indices)
        : m_front(valid_indices.front()...)
        , m_size(valid_indices.extents()...)
        , m_stride(valid_indices.strides()...)
    {
    }

    /**
     * @brief Check if the specified DiscreteElement is found in this subdomain.
     *
     * @param elem The element which may or may not be in this subdomain.
     *
     * @returns bool True if the element is found in this subdomain, False otherwise.
     */
    template <class... DDims>
    KOKKOS_FUNCTION bool contains(ddc::DiscreteElement<DDims...> elem) const
    {
        static_assert(
                (ddc::in_tags_v<DDims, tags_seq> && ...),
                "Requested Tag absent from DiscreteSubDomain");

        ddc::DiscreteVector<DDims...> distance(elem - ddc::select<DDims...>(m_front));
        return (((ddc::get<DDims>(distance) % ddc::get<DDims>(m_stride)) == 0) && ...);
    }

    /**
     * @brief Check if all elements of the specified DiscreteDomain are found in this subdomain.
     *
     * @param dom The domain which may or may not be in this subdomain.
     *
     * @returns bool True if the all elements of the specified DiscreteDomain are found in this subdomain.
     *               False otherwise.
     */
    template <class... DDims>
    KOKKOS_FUNCTION bool contains(ddc::DiscreteDomain<DDims...> dom) const
    {
        static_assert(
                (ddc::in_tags_v<DDims, tags_seq> && ...),
                "Requested Tag absent from DiscreteSubDomain");
        return contains(dom.front())
               && (((dom.template extent<DDims>().value() == 1)
                    || (ddc::select<DDims>(m_stride) == 1))
                   && ...);
    }

    /**
     * @brief Get the index of the DiscreteElement within the subdomain.
     * This function is particularly useful to index an mdspan over the subdomain.
     *
     * @param elem A 1D DiscreteElement which is inside the subdomain.
     *
     * @returns The index of the element.
     */
    template <class Dim, class = std::enable_if_t<ddc::in_tags_v<Dim, tags_seq>>>
    KOKKOS_FUNCTION std::size_t get_index(ddc::DiscreteElement<Dim> elem) const
    {
        assert(contains(elem));
        return (elem - ddc::select<Dim>(m_front)).value() / ddc::get<Dim>(m_stride);
    }

    /**
     * @brief Get the size of the subdomain in each dimension.
     *
     * @returns A DiscreteVector describing the size of the subdomain in each dimension.
     */
    KOKKOS_FUNCTION constexpr ddc::DiscreteVector<Dims...> extents() const noexcept
    {
        return m_size;
    }

    /**
     * @brief Get the size of the subdomain in the specified dimension.
     *
     * @tparam QueryDDim The dimension in which the size is requested.
     *
     * @returns A DiscreteVector describing the size of the subdomain in the specified dimension.
     */
    template <class QueryDDim>
    KOKKOS_FUNCTION constexpr ddc::DiscreteVector<QueryDDim> extent() const noexcept
    {
        return ddc::select<QueryDDim>(m_size);
    }

    /**
     * @brief Get the total number of elements in the subdomain.
     *
     * @return The total number of elements in the subdomain.
     */
    KOKKOS_FUNCTION constexpr std::size_t size() const
    {
        return (1ul * ... * (ddc::get<Dims>(m_size)));
    }

    /**
     * @brief Get the first element in the subdomain.
     *
     * @return The first element in the subdomain.
     */
    KOKKOS_FUNCTION constexpr ddc::DiscreteElement<Dims...> front() const noexcept
    {
        return m_front;
    }

    /**
     * @brief Get the last element in the subdomain.
     *
     * @return The last element in the subdomain.
     */
    KOKKOS_FUNCTION constexpr ddc::DiscreteElement<Dims...> back() const noexcept
    {
        return m_front + m_stride * (m_size - 1);
    }

    /**
     * @brief Get the stride from one element of the subdomain to another.
     *
     * @tparam QueryDim The dimension being queried
     *
     * @return A 1D DiscreteVector describing a stride from one element to another.
     */
    template <class QueryDim>
    KOKKOS_FUNCTION constexpr ddc::DiscreteVector<QueryDim> stride() const noexcept
    {
        return ddc::select<QueryDim>(m_stride);
    }

    /**
     * @brief Get the strides from one element of the subdomain to another.
     *
     * @return A DiscreteVector describing the strides from one element to another.
     */
    KOKKOS_FUNCTION constexpr ddc::DiscreteVector<Dims...> strides() const noexcept
    {
        return m_stride;
    }

    /**
     * @brief Get the iterator to the first element of the DiscreteSubDomain.
     * The elements are pairs of DiscreteElements and indices.
     *
     * @returns The iterator to the first element of the DiscreteSubDomain.
     */
    KOKKOS_FUNCTION auto begin() const
    {
        static_assert(sizeof...(Dims) == 1);
        return DiscreteSubDomainIterator(m_front, m_stride);
    }

    /**
     * @brief Get the iterator to the end of the DiscreteSubDomain.
     * The elements are pairs of DiscreteElements and indices.
     *
     * @returns The iterator to the end of the DiscreteSubDomain.
     */
    KOKKOS_FUNCTION auto end() const
    {
        static_assert(sizeof...(Dims) == 1);
        return DiscreteSubDomainIterator(m_front + m_stride * m_size, m_stride);
    }
};

/**
 * @brief An iterator type for the DiscreteSubDomain.
 */
template <class Grid1D>
struct DiscreteSubDomainIterator
{
private:
    ddc::DiscreteElement<Grid1D> m_value;
    ddc::DiscreteVector<Grid1D> m_stride;

public:
    /// The type of iterator
    using iterator_category = std::random_access_iterator_tag;

    /// The type of the values stored in the iterator
    using value_type = ddc::DiscreteElement<Grid1D>;

    /// The type of the stride between values
    using stride_type = ddc::DiscreteVector<Grid1D>;

    /// The type that can be used to increment the iterator
    using difference_type = std::ptrdiff_t;

    KOKKOS_DEFAULTED_FUNCTION DiscreteSubDomainIterator() = default;

    /**
     * @brief Build an iterator using the current value and the distance to the following element
     *
     * @param value The value of the discrete sub-domain element.
     * @param stride The stride between consectuive sub-domain elements.
     */
    KOKKOS_FUNCTION constexpr explicit DiscreteSubDomainIterator(
            ddc::DiscreteElement<Grid1D> value,
            ddc::DiscreteVector<Grid1D> stride)
        : m_value(value)
        , m_stride(stride)
    {
    }

    /**
     * @brief Get the value referred to by the iterator.
     * @return The value referred to by the iterator.
     */
    KOKKOS_FUNCTION constexpr ddc::DiscreteElement<Grid1D> operator*() const noexcept
    {
        return m_value;
    }

    /**
     * @brief The prefix increment operator.
     * @return A reference to the current incremented iterator.
     */
    KOKKOS_FUNCTION constexpr DiscreteSubDomainIterator& operator++()
    {
        m_value = m_value + m_stride;
        return *this;
    }

    /**
     * @brief The postfix increment operator.
     * @return A reference to the non-incremented iterator.
     */
    KOKKOS_FUNCTION constexpr DiscreteSubDomainIterator operator++(int)
    {
        auto tmp = *this;
        ++*this;
        return tmp;
    }

    /**
     * @brief The prefix decrement operator.
     * @return A reference to the current decremented iterator.
     */
    KOKKOS_FUNCTION constexpr DiscreteSubDomainIterator& operator--()
    {
        m_value = m_value - m_stride;
        return *this;
    }

    /**
     * @brief The postfix decrement operator.
     * @return A reference to the non-decremented iterator.
     */
    KOKKOS_FUNCTION constexpr DiscreteSubDomainIterator operator--(int)
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
    KOKKOS_FUNCTION constexpr DiscreteSubDomainIterator& operator+=(difference_type n)
    {
        m_value = m_value + n * m_stride;
        return *this;
    }

    /**
     * @brief Decrement the current iterator by n elements.
     * @param n The number of strides by which the iterator should be decremented.
     * @return A reference to the current decremented iterator.
     */
    KOKKOS_FUNCTION constexpr DiscreteSubDomainIterator& operator-=(difference_type n)
    {
        m_value = m_value - n * m_stride;
        return *this;
    }

    /**
     * @brief Access the n-th following element.
     * @param n The index of the element.
     * @return The DiscreteElement n-th position in the sub-domain following the value
     *          indicated by this iterator.
     */
    KOKKOS_FUNCTION constexpr ddc::DiscreteElement<Grid1D> operator[](difference_type n) const
    {
        return m_value + n * m_stride;
    }

    /// Compare iterator equality
    friend KOKKOS_FUNCTION constexpr bool operator==(
            DiscreteSubDomainIterator const& xx,
            DiscreteSubDomainIterator const& yy)
    {
        return xx.m_value == yy.m_value;
    }

    /// Compare iterator non-equality
    friend KOKKOS_FUNCTION constexpr bool operator!=(
            DiscreteSubDomainIterator const& xx,
            DiscreteSubDomainIterator const& yy)
    {
        return xx.m_value != yy.m_value;
    }

    /// Compare the order of iterators
    friend KOKKOS_FUNCTION constexpr bool operator<(
            DiscreteSubDomainIterator const& xx,
            DiscreteSubDomainIterator const& yy)
    {
        return xx.m_value < yy.m_value;
    }

    /// Compare the order of iterators
    friend KOKKOS_FUNCTION constexpr bool operator>(
            DiscreteSubDomainIterator const& xx,
            DiscreteSubDomainIterator const& yy)
    {
        return yy < xx;
    }

    /// Compare the order of iterators
    friend KOKKOS_FUNCTION constexpr bool operator<=(
            DiscreteSubDomainIterator const& xx,
            DiscreteSubDomainIterator const& yy)
    {
        return !(yy < xx);
    }

    /// Compare the order of iterators
    friend KOKKOS_FUNCTION constexpr bool operator>=(
            DiscreteSubDomainIterator const& xx,
            DiscreteSubDomainIterator const& yy)
    {
        return !(xx < yy);
    }

    /// Increment an iterator by n elements.
    friend KOKKOS_FUNCTION constexpr DiscreteSubDomainIterator operator+(
            DiscreteSubDomainIterator i,
            difference_type n)
    {
        return i += n;
    }

    /// Increment an iterator by n elements.
    friend KOKKOS_FUNCTION constexpr DiscreteSubDomainIterator operator+(
            difference_type n,
            DiscreteSubDomainIterator i)
    {
        return i += n;
    }

    /// Decrement an iterator by n elements.
    friend KOKKOS_FUNCTION constexpr DiscreteSubDomainIterator operator-(
            DiscreteSubDomainIterator i,
            difference_type n)
    {
        return i -= n;
    }

    /// Decrement an iterator by n elements.
    friend KOKKOS_FUNCTION constexpr difference_type operator-(
            DiscreteSubDomainIterator const& xx,
            DiscreteSubDomainIterator const& yy)
    {
        assert(xx.stride() == yy.stride());
        return (yy.m_value > xx.m_value)
                       ? ((-static_cast<difference_type>(yy.m_value - xx.m_value))
                          / ddc::get<Grid1D>(xx.m_stride))
                       : ((xx.m_value - yy.m_value) / ddc::get<Grid1D>(xx.m_stride));
    }
};

/// A class to create a DiscreteSubDomain type from a TypeSeq.
template <class T>
struct ToDiscreteSubDomain;

template <class... Dims>
struct ToDiscreteSubDomain<ddc::detail::TypeSeq<Dims...>>
{
    using value = DiscreteSubDomain<Dims...>;
};

/// @brief A templated type for creating a DiscreteSubDomain from a TypeSeq.
template <class T>
using to_subdomain_collection = typename ToDiscreteSubDomain<T>::value;
