// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <tuple>

#include "ddc/coordinate.hpp"
#include "ddc/detail/type_seq.hpp"
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_vector.hpp"

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
    using mcoord_type = DiscreteCoordinate<DDims...>;

    using mlength_type = DiscreteVector<DDims...>;

private:
    DiscreteCoordinate<DDims...> m_lbound;

    DiscreteCoordinate<DDims...> m_ubound;

public:
    static constexpr std::size_t rank()
    {
        return (0 + ... + DDims::rank());
    }

    DiscreteDomain() = default;

    /// Construct a DiscreteDomain from a reordered copy of `domain`
    template <class... ODDims>
    explicit constexpr DiscreteDomain(DiscreteDomain<ODDims...> const& domain)
        : m_lbound(domain.front())
        , m_ubound(domain.back())
    {
    }

    // Use SFINAE to disambiguate with the copy constructor.
    // Note that SFINAE may be redundant because a template constructor should not be selected as a copy constructor.
    template <std::size_t N = sizeof...(DDims), class = std::enable_if_t<(N != 1)>>
    explicit constexpr DiscreteDomain(DiscreteDomain<DDims> const&... domains)
        : m_lbound(domains.front()...)
        , m_ubound(domains.back()...)
    {
    }

    /** Construct a DiscreteDomain starting from (0, ..., 0) with size points.
     * @param size the number of points in each dimension
     * 
     * @deprecated use the version with explicit lower bound instead
     */
    [[deprecated]] constexpr DiscreteDomain(mlength_type const& size)
        : m_lbound((get<DDims>(size) - get<DDims>(size))...) // Hack to have expansion of zero
        , m_ubound((get<DDims>(size) - 1)...)
    {
    }

    /** Construct a DiscreteDomain starting from lbound with size points.
     * @param lbound the lower bound in each direction
     * @param size the number of points in each direction
     */
    constexpr DiscreteDomain(mcoord_type const& lbound, mlength_type const& size)
        : m_lbound(lbound)
        , m_ubound((uid<DDims>(lbound) + get<DDims>(size) - 1)...)
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
        return m_lbound == other.m_lbound && m_ubound == other.m_ubound;
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
        return (1ul * ... * (uid<DDims>(m_ubound) + 1 - uid<DDims>(m_lbound)));
    }

    constexpr mlength_type extents() const noexcept
    {
        return mlength_type((uid<DDims>(m_ubound) + 1 - uid<DDims>(m_lbound))...);
    }

    template <class QueryDDim>
    inline constexpr DiscreteVector<QueryDDim> extent() const noexcept
    {
        return DiscreteVector<QueryDDim>(uid<QueryDDim>(m_ubound) + 1 - uid<QueryDDim>(m_lbound));
    }

    constexpr mcoord_type front() const noexcept
    {
        return m_lbound;
    }

    constexpr mcoord_type back() const noexcept
    {
        return m_ubound;
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
        assert(((uid<ODDims>(m_lbound) <= uid<ODDims>(odomain.m_lbound)) && ...));
        assert(((uid<ODDims>(m_ubound) >= uid<ODDims>(odomain.m_ubound)) && ...));
        const DiscreteVector<DDims...> myextents = extents();
        const DiscreteVector<ODDims...> oextents = odomain.extents();
        return DiscreteDomain(
                DiscreteCoordinate<DDims...>(
                        (uid_or<DDims>(odomain.m_lbound, uid<DDims>(m_lbound)))...),
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
        return DiscreteDomainIterator<DDim0>(DiscreteCoordinate<DDims...>(back() + 1));
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
        return DiscreteDomainIterator<DDim0>(DiscreteCoordinate<DDims...>(back() + 1));
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
constexpr DiscreteCoordinate<QueryDDims...> front(DiscreteDomain<DDims...> const& domain) noexcept
{
    return DiscreteCoordinate<QueryDDims...>(select<QueryDDims>(domain).front()...);
}

template <class... QueryDDims, class... DDims>
constexpr DiscreteCoordinate<QueryDDims...> back(DiscreteDomain<DDims...> const& domain) noexcept
{
    return DiscreteCoordinate<QueryDDims...>(select<QueryDDims>(domain).back()...);
}

template <class... QueryDDims, class... DDims>
Coordinate<QueryDDims...> to_real(
        DiscreteDomain<DDims...> const& domain,
        DiscreteCoordinate<QueryDDims...> const& icoord) noexcept
{
    return Coordinate<QueryDDims...>(
            select<QueryDDims>(domain).to_real(select<QueryDDims>(icoord))...);
}

template <class... QueryDDims, class... DDims>
Coordinate<QueryDDims...> rmin(DiscreteDomain<DDims...> const& domain) noexcept
{
    return Coordinate<QueryDDims...>(select<QueryDDims>(domain).rmin()...);
}

template <class... QueryDDims, class... DDims>
Coordinate<QueryDDims...> rmax(DiscreteDomain<DDims...> const& domain) noexcept
{
    return Coordinate<QueryDDims...>(select<QueryDDims>(domain).rmax()...);
}

namespace detail {

template <class QueryDDimSeq>
struct Selection;

template <class... QueryDDims>
struct Selection<detail::TypeSeq<QueryDDims...>>
{
    template <class Domain>
    static constexpr auto select(Domain const& domain)
    {
        return ::select<QueryDDims...>(domain);
    }
};

} // namespace detail

template <class QueryDDimSeq, class... DDims>
constexpr auto select_by_type_seq(DiscreteDomain<DDims...> const& domain)
{
    return detail::Selection<QueryDDimSeq>::select(domain);
}

template <class DDim>
struct DiscreteDomainIterator
{
private:
    DiscreteCoordinate<DDim> m_value = DiscreteCoordinate<DDim>();

public:
    using iterator_category = std::random_access_iterator_tag;

    using value_type = DiscreteCoordinate<DDim>;

    using difference_type = std::ptrdiff_t;

    DiscreteDomainIterator() = default;

    constexpr explicit DiscreteDomainIterator(DiscreteCoordinate<DDim> __value) : m_value(__value)
    {
    }

    constexpr DiscreteCoordinate<DDim> operator*() const noexcept
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
            m_value.uid() += static_cast<DiscreteCoordinateElement>(__n);
        else
            m_value.uid() -= static_cast<DiscreteCoordinateElement>(-__n);
        return *this;
    }

    constexpr DiscreteDomainIterator& operator-=(difference_type __n)
    {
        if (__n >= difference_type(0))
            m_value.uid() -= static_cast<DiscreteCoordinateElement>(__n);
        else
            m_value.uid() += static_cast<DiscreteCoordinateElement>(-__n);
        return *this;
    }

    constexpr DiscreteCoordinate<DDim> operator[](difference_type __n) const
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
