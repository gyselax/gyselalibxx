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
        return sizeof...(DDims);
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
        return DiscreteDomain(front() + (extents() - n), n);
    }

    constexpr DiscreteDomain remove_first(mlength_type n) const
    {
        return DiscreteDomain(front() + n, extents() - n);
    }

    constexpr DiscreteDomain remove_last(mlength_type n) const
    {
        return DiscreteDomain(front(), extents() - n);
    }

    constexpr DiscreteDomain remove(mlength_type n1, mlength_type n2) const
    {
        return DiscreteDomain(front() + n1, extents() - n1 - n2);
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
    constexpr decltype(auto) operator[](std::size_t n)
    {
        return begin()[n];
    }

    template <
            std::size_t N = sizeof...(DDims),
            class = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<DDims...>>>>
    constexpr decltype(auto) operator[](std::size_t n) const
    {
        return begin()[n];
    }
};

template <>
class DiscreteDomain<>
{
    template <class...>
    friend class DiscreteDomain;

public:
    using discrete_element_type = DiscreteElement<>;

    using mlength_type = DiscreteVector<>;

    static constexpr std::size_t rank()
    {
        return 0;
    }

    constexpr DiscreteDomain() = default;

    // Construct a DiscreteDomain from a reordered copy of `domain`
    template <class... ODDims>
    explicit constexpr DiscreteDomain([[maybe_unused]] DiscreteDomain<ODDims...> const& domain)
    {
    }

    /** Construct a DiscreteDomain starting from element_begin with size points.
     * @param element_begin the lower bound in each direction
     * @param size the number of points in each direction
     */
    constexpr DiscreteDomain(
            [[maybe_unused]] discrete_element_type const& element_begin,
            [[maybe_unused]] mlength_type const& size)
    {
    }

    constexpr DiscreteDomain(DiscreteDomain const& x) = default;

    constexpr DiscreteDomain(DiscreteDomain&& x) = default;

    ~DiscreteDomain() = default;

    DiscreteDomain& operator=(DiscreteDomain const& x) = default;

    DiscreteDomain& operator=(DiscreteDomain&& x) = default;

    constexpr bool operator==([[maybe_unused]] DiscreteDomain const& other) const
    {
        return true;
    }

#if __cplusplus <= 201703L
    // Shall not be necessary anymore in C++20
    // `a!=b` shall be translated by the compiler to `!(a==b)`
    constexpr bool operator!=(DiscreteDomain const& other) const
    {
        return !(*this == other);
    }
#endif

    constexpr std::size_t size() const
    {
        return 1;
    }

    constexpr mlength_type extents() const noexcept
    {
        return {};
    }

    constexpr discrete_element_type front() const noexcept
    {
        return {};
    }

    constexpr discrete_element_type back() const noexcept
    {
        return {};
    }

    constexpr DiscreteDomain take_first([[maybe_unused]] mlength_type n) const
    {
        return *this;
    }

    constexpr DiscreteDomain take_last([[maybe_unused]] mlength_type n) const
    {
        return *this;
    }

    constexpr DiscreteDomain remove_first([[maybe_unused]] mlength_type n) const
    {
        return *this;
    }

    constexpr DiscreteDomain remove_last([[maybe_unused]] mlength_type n) const
    {
        return *this;
    }

    constexpr DiscreteDomain remove(
            [[maybe_unused]] mlength_type n1,
            [[maybe_unused]] mlength_type n2) const
    {
        return *this;
    }

    template <class... ODims>
    constexpr DiscreteDomain restrict(DiscreteDomain<ODims...> const&) const
    {
        return *this;
    }

    constexpr bool empty() const noexcept
    {
        return false;
    }

    constexpr explicit operator bool()
    {
        return true;
    }
};

template <class... QueryDDims, class... DDims>
constexpr DiscreteDomain<QueryDDims...> select(DiscreteDomain<DDims...> const& domain)
{
    return DiscreteDomain<QueryDDims...>(
            select<QueryDDims...>(domain.front()),
            select<QueryDDims...>(domain.extents()));
}

namespace detail {

template <class T>
struct ConvertTypeSeqToDiscreteDomain;

template <class... DDims>
struct ConvertTypeSeqToDiscreteDomain<detail::TypeSeq<DDims...>>
{
    using type = DiscreteDomain<DDims...>;
};

template <class T>
using convert_type_seq_to_discrete_domain = typename ConvertTypeSeqToDiscreteDomain<T>::type;

} // namespace detail

// Computes the cartesian product of DiscreteDomain types
// Example usage : "using DDom = cartesian_prod_t<DDom1,DDom2,DDom3>;"
template <typename... DDom>
struct cartesian_prod;

template <typename... DDim1, typename... DDim2>
struct cartesian_prod<ddc::DiscreteDomain<DDim1...>, ddc::DiscreteDomain<DDim2...>>
{
    using type = ddc::DiscreteDomain<DDim1..., DDim2...>;
};

template <typename DDom>
struct cartesian_prod<DDom>
{
    using type = DDom;
};

template <typename DDom1, typename DDom2, typename... Tail>
struct cartesian_prod<DDom1, DDom2, Tail...>
{
    using type =
            typename cartesian_prod<typename cartesian_prod<DDom1, DDom2>::type, Tail...>::type;
};

template <typename... DDom>
using cartesian_prod_t = typename cartesian_prod<DDom...>::type;

// Computes the substraction DDom_a - DDom_b in the sense of linear spaces(retained dimensions are those in DDom_a which are not in DDom_b)
template <class... DDimsA, class... DDimsB>
constexpr auto remove_dims_of(
        DiscreteDomain<DDimsA...> const& DDom_a,
        [[maybe_unused]] DiscreteDomain<DDimsB...> const& DDom_b) noexcept
{
    using TagSeqA = detail::TypeSeq<DDimsA...>;
    using TagSeqB = detail::TypeSeq<DDimsB...>;

    using type_seq_r = type_seq_remove_t<TagSeqA, TagSeqB>;
    return detail::convert_type_seq_to_discrete_domain<type_seq_r>(DDom_a);
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

namespace detail {

template <class QueryDDimSeq>
struct Selection;

template <class... QueryDDims>
struct Selection<detail::TypeSeq<QueryDDims...>>
{
    template <class Domain>
    static constexpr auto select(Domain const& domain)
    {
        return ddc::select<QueryDDims...>(domain);
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
    DiscreteElement<DDim> m_value = DiscreteElement<DDim>();

public:
    using iterator_category = std::random_access_iterator_tag;

    using value_type = DiscreteElement<DDim>;

    using difference_type = std::ptrdiff_t;

    DiscreteDomainIterator() = default;

    constexpr explicit DiscreteDomainIterator(DiscreteElement<DDim> value) : m_value(value) {}

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
        auto tmp = *this;
        ++*this;
        return tmp;
    }

    constexpr DiscreteDomainIterator& operator--()
    {
        --m_value.uid();
        return *this;
    }

    constexpr DiscreteDomainIterator operator--(int)
    {
        auto tmp = *this;
        --*this;
        return tmp;
    }

    constexpr DiscreteDomainIterator& operator+=(difference_type n)
    {
        if (n >= difference_type(0))
            m_value.uid() += static_cast<DiscreteElementType>(n);
        else
            m_value.uid() -= static_cast<DiscreteElementType>(-n);
        return *this;
    }

    constexpr DiscreteDomainIterator& operator-=(difference_type n)
    {
        if (n >= difference_type(0))
            m_value.uid() -= static_cast<DiscreteElementType>(n);
        else
            m_value.uid() += static_cast<DiscreteElementType>(-n);
        return *this;
    }

    constexpr DiscreteElement<DDim> operator[](difference_type n) const
    {
        return m_value + n;
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

    friend constexpr DiscreteDomainIterator operator+(DiscreteDomainIterator i, difference_type n)
    {
        return i += n;
    }

    friend constexpr DiscreteDomainIterator operator+(difference_type n, DiscreteDomainIterator i)
    {
        return i += n;
    }

    friend constexpr DiscreteDomainIterator operator-(DiscreteDomainIterator i, difference_type n)
    {
        return i -= n;
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
