// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <tuple>

#include "ddc/coordinate.hpp"
#include "ddc/detail/type_seq.hpp"
#include "ddc/discrete_element.hpp"
#include "ddc/discrete_vector.hpp"

namespace ddc {

template <class DDim>
struct DiscreteDomainIterator;

template <class... DDims>
class DiscreteDomain;

template <class... DDims>
class DiscreteDomain
{
    template <class...>
    friend class DiscreteDomain;

public:
    using discrete_element_type = DiscreteElement<DDims...>;

    using mlength_type = DiscreteVector<DDims...>;

private:
    DiscreteElement<DDims...> m_element_begin;

    DiscreteElement<DDims...> m_element_end;

public:
    static constexpr std::size_t rank()
    {
        return (0 + ... + DDims::rank());
    }

    DiscreteDomain() = default;

    /// Construct a DiscreteDomain from a reordered copy of `domain`
    template <class... ODDims>
    explicit constexpr DiscreteDomain(DiscreteDomain<ODDims...> const& domain)
        : m_element_begin(domain.front())
        , m_element_end(domain.front() + domain.extents())
    {
    }

    // Use SFINAE to disambiguate with the copy constructor.
    // Note that SFINAE may be redundant because a template constructor should not be selected as a copy constructor.
    template <std::size_t N = sizeof...(DDims), class = std::enable_if_t<(N != 1)>>
    explicit constexpr DiscreteDomain(DiscreteDomain<DDims> const&... domains)
        : m_element_begin(domains.front()...)
        , m_element_end((domains.front() + domains.extents())...)
    {
    }

    /** Construct a DiscreteDomain starting from (0, ..., 0) with size points.
     * @param size the number of points in each dimension
     * 
     * @deprecated use the version with explicit lower bound instead
     */
    [[deprecated]] constexpr DiscreteDomain(mlength_type const& size)
        : m_element_begin(
                (get<DDims>(size) - get<DDims>(size))...) // Hack to have expansion of zero
        , m_element_end(get<DDims>(size)...)
    {
    }

    /** Construct a DiscreteDomain starting from element_begin with size points.
     * @param element_begin the lower bound in each direction
     * @param size the number of points in each direction
     */
    constexpr DiscreteDomain(discrete_element_type const& element_begin, mlength_type const& size)
        : m_element_begin(element_begin)
        , m_element_end(element_begin + size)
    {
    }

    DiscreteDomain(DiscreteDomain const& x) = default;

    DiscreteDomain(DiscreteDomain&& x) = default;

    ~DiscreteDomain() = default;

    DiscreteDomain& operator=(DiscreteDomain const& x) = default;

    DiscreteDomain& operator=(DiscreteDomain&& x) = default;

    template <class... ODims>
    constexpr bool operator==(DiscreteDomain<ODims...> const& other) const
    {
        return m_element_begin == other.m_element_begin && m_element_end == other.m_element_end;
    }

#if __cplusplus <= 201703L
    // Shall not be necessary anymore in C++20
    // `a!=b` shall be translated by the compiler to `!(a==b)`
    constexpr bool operator!=(DiscreteDomain const& other) const
    {
        return !(*this == other);
    }
#endif

    std::size_t size() const
    {
        return (1ul * ... * (uid<DDims>(m_element_end) - uid<DDims>(m_element_begin)));
    }

    constexpr mlength_type extents() const noexcept
    {
        return mlength_type((uid<DDims>(m_element_end) - uid<DDims>(m_element_begin))...);
    }

    template <class QueryDDim>
    inline constexpr DiscreteVector<QueryDDim> extent() const noexcept
    {
        return DiscreteVector<QueryDDim>(
                uid<QueryDDim>(m_element_end) - uid<QueryDDim>(m_element_begin));
    }

    constexpr discrete_element_type front() const noexcept
    {
        return m_element_begin;
    }

    constexpr discrete_element_type back() const noexcept
    {
        return discrete_element_type((uid<DDims>(m_element_end) - 1)...);
    }

    constexpr DiscreteDomain take_first(mlength_type n) const
    {
        return DiscreteDomain(front(), n);
    }

    constexpr DiscreteDomain take_last(mlength_type n) const
    {
        return DiscreteDomain(front() + size() - n, n);
    }

    constexpr DiscreteDomain remove_first(mlength_type n) const
    {
        return take_last(size() - n);
    }

    constexpr DiscreteDomain remove_last(mlength_type n) const
    {
        return take_first(size() - n);
    }

    constexpr DiscreteDomain remove(mlength_type n1, mlength_type n2) const
    {
        return remove_first(n1).remove_last(n2);
    }

    template <class... ODDims>
    constexpr auto restrict(DiscreteDomain<ODDims...> const& odomain) const
    {
        assert(((uid<ODDims>(m_element_begin) <= uid<ODDims>(odomain.m_element_begin)) && ...));
        assert(((uid<ODDims>(m_element_end) >= uid<ODDims>(odomain.m_element_end)) && ...));
        const DiscreteVector<DDims...> myextents = extents();
        const DiscreteVector<ODDims...> oextents = odomain.extents();
        return DiscreteDomain(
                DiscreteElement<DDims...>(
                        (uid_or<DDims>(odomain.m_element_begin, uid<DDims>(m_element_begin)))...),
                DiscreteVector<DDims...>((get_or<DDims>(oextents, get<DDims>(myextents)))...));
    }

    constexpr bool empty() const noexcept
    {
        return size() == 0;
    }

    constexpr explicit operator bool()
    {
        return !empty();
    }

    template <
            std::size_t N = sizeof...(DDims),
            class DDim0 = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<DDims...>>>>
    auto begin() const
    {
        return DiscreteDomainIterator<DDim0>(front());
    }

    template <
            std::size_t N = sizeof...(DDims),
            class DDim0 = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<DDims...>>>>
    auto end() const
    {
        return DiscreteDomainIterator<DDim0>(m_element_end);
    }

    template <
            std::size_t N = sizeof...(DDims),
            class DDim0 = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<DDims...>>>>
    auto cbegin() const
    {
        return DiscreteDomainIterator<DDim0>(front());
    }

    template <
            std::size_t N = sizeof...(DDims),
            class DDim0 = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<DDims...>>>>
    auto cend() const
    {
        return DiscreteDomainIterator<DDim0>(m_element_end);
    }

    template <
            std::size_t N = sizeof...(DDims),
            class = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<DDims...>>>>
    constexpr decltype(auto) operator[](std::size_t __n)
    {
        return begin()[__n];
    }

    template <
            std::size_t N = sizeof...(DDims),
            class = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<DDims...>>>>
    constexpr decltype(auto) operator[](std::size_t __n) const
    {
        return begin()[__n];
    }
};

template <class... QueryDDims, class... DDims>
constexpr DiscreteDomain<QueryDDims...> select(DiscreteDomain<DDims...> const& domain)
{
    return DiscreteDomain<QueryDDims...>(
            select<QueryDDims...>(domain.front()),
            select<QueryDDims...>(domain.extents()));
}

template <class... QueryDDims, class... DDims>
constexpr DiscreteVector<QueryDDims...> extents(DiscreteDomain<DDims...> const& domain) noexcept
{
    return DiscreteVector<QueryDDims...>(select<QueryDDims>(domain).size()...);
}

template <class... QueryDDims, class... DDims>
constexpr DiscreteElement<QueryDDims...> front(DiscreteDomain<DDims...> const& domain) noexcept
{
    return DiscreteElement<QueryDDims...>(select<QueryDDims>(domain).front()...);
}

template <class... QueryDDims, class... DDims>
constexpr DiscreteElement<QueryDDims...> back(DiscreteDomain<DDims...> const& domain) noexcept
{
    return DiscreteElement<QueryDDims...>(select<QueryDDims>(domain).back()...);
}

template <class... QueryDDims, class... DDims>
ddc::Coordinate<QueryDDims...> coordinate(
        DiscreteDomain<DDims...> const& domain,
        DiscreteElement<QueryDDims...> const& icoord) noexcept
{
    return ddc::Coordinate<QueryDDims...>(
            select<QueryDDims>(domain).coordinate(select<QueryDDims>(icoord))...);
}

template <class... QueryDDims, class... DDims>
ddc::Coordinate<QueryDDims...> rmin(DiscreteDomain<DDims...> const& domain) noexcept
{
    return ddc::Coordinate<QueryDDims...>(select<QueryDDims>(domain).rmin()...);
}

template <class... QueryDDims, class... DDims>
ddc::Coordinate<QueryDDims...> rmax(DiscreteDomain<DDims...> const& domain) noexcept
{
    return ddc::Coordinate<QueryDDims...>(select<QueryDDims>(domain).rmax()...);
}

namespace ddc_detail {

template <class QueryDDimSeq>
struct Selection;

template <class... QueryDDims>
struct Selection<ddc_detail::TypeSeq<QueryDDims...>>
{
    template <class Domain>
    static constexpr auto select(Domain const& domain)
    {
        return ddc::select<QueryDDims...>(domain);
    }
};

} // namespace ddc_detail

template <class QueryDDimSeq, class... DDims>
constexpr auto select_by_type_seq(DiscreteDomain<DDims...> const& domain)
{
    return ddc_detail::Selection<QueryDDimSeq>::select(domain);
}

template <class DDim>
struct DiscreteDomainIterator
{
private:
    DiscreteElement<DDim> m_value = DiscreteElement<DDim>();

public:
    using iterator_category = std::random_access_iterator_tag;

    using value_type = DiscreteElement<DDim>;

    using difference_type = std::ptrdiff_t;

    DiscreteDomainIterator() = default;

    constexpr explicit DiscreteDomainIterator(DiscreteElement<DDim> __value) : m_value(__value) {}

    constexpr DiscreteElement<DDim> operator*() const noexcept
    {
        return m_value;
    }

    constexpr DiscreteDomainIterator& operator++()
    {
        ++m_value.uid();
        return *this;
    }

    constexpr DiscreteDomainIterator operator++(int)
    {
        auto __tmp = *this;
        ++*this;
        return __tmp;
    }

    constexpr DiscreteDomainIterator& operator--()
    {
        --m_value.uid();
        return *this;
    }

    constexpr DiscreteDomainIterator operator--(int)
    {
        auto __tmp = *this;
        --*this;
        return __tmp;
    }

    constexpr DiscreteDomainIterator& operator+=(difference_type __n)
    {
        if (__n >= difference_type(0))
            m_value.uid() += static_cast<DiscreteElementType>(__n);
        else
            m_value.uid() -= static_cast<DiscreteElementType>(-__n);
        return *this;
    }

    constexpr DiscreteDomainIterator& operator-=(difference_type __n)
    {
        if (__n >= difference_type(0))
            m_value.uid() -= static_cast<DiscreteElementType>(__n);
        else
            m_value.uid() += static_cast<DiscreteElementType>(-__n);
        return *this;
    }

    constexpr DiscreteElement<DDim> operator[](difference_type __n) const
    {
        return m_value + __n;
    }

    friend constexpr bool operator==(
            DiscreteDomainIterator const& xx,
            DiscreteDomainIterator const& yy)
    {
        return xx.m_value == yy.m_value;
    }

    friend constexpr bool operator!=(
            DiscreteDomainIterator const& xx,
            DiscreteDomainIterator const& yy)
    {
        return xx.m_value != yy.m_value;
    }

    friend constexpr bool operator<(
            DiscreteDomainIterator const& xx,
            DiscreteDomainIterator const& yy)
    {
        return xx.m_value < yy.m_value;
    }

    friend constexpr bool operator>(
            DiscreteDomainIterator const& xx,
            DiscreteDomainIterator const& yy)
    {
        return yy < xx;
    }

    friend constexpr bool operator<=(
            DiscreteDomainIterator const& xx,
            DiscreteDomainIterator const& yy)
    {
        return !(yy < xx);
    }

    friend constexpr bool operator>=(
            DiscreteDomainIterator const& xx,
            DiscreteDomainIterator const& yy)
    {
        return !(xx < yy);
    }

    friend constexpr DiscreteDomainIterator operator+(
            DiscreteDomainIterator __i,
            difference_type __n)
    {
        return __i += __n;
    }

    friend constexpr DiscreteDomainIterator operator+(
            difference_type __n,
            DiscreteDomainIterator __i)
    {
        return __i += __n;
    }

    friend constexpr DiscreteDomainIterator operator-(
            DiscreteDomainIterator __i,
            difference_type __n)
    {
        return __i -= __n;
    }

    friend constexpr difference_type operator-(
            DiscreteDomainIterator const& xx,
            DiscreteDomainIterator const& yy)
    {
        return (yy.m_value > xx.m_value) ? (-static_cast<difference_type>(yy.m_value - xx.m_value))
                                         : (xx.m_value - yy.m_value);
    }
};

} // namespace ddc
