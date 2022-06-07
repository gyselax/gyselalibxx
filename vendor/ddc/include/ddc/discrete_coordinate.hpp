// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstddef>
#include <ostream>
#include <utility>

#include "ddc/detail/type_seq.hpp"
#include "ddc/discrete_vector.hpp"


template <class...>
class DiscreteCoordinate;

template <class T>
struct IsDiscreteCoordinate : std::false_type
{
};

template <class... Tags>
struct IsDiscreteCoordinate<DiscreteCoordinate<Tags...>> : std::true_type
{
};

template <class T>
inline constexpr bool is_discrete_coordinate_v = IsDiscreteCoordinate<T>::value;


/** A DiscreteCoordElement is a scalar that identifies an element of the discrete dimension
 */
using DiscreteCoordinateElement = std::size_t;

template <class Tag>
inline constexpr DiscreteCoordinateElement const& uid(DiscreteCoordinate<Tag> const& tuple) noexcept
{
    return tuple.uid();
}

template <class Tag>
inline constexpr DiscreteCoordinateElement& uid(DiscreteCoordinate<Tag>& tuple) noexcept
{
    return tuple.uid();
}

template <class QueryTag, class... Tags>
inline constexpr DiscreteCoordinateElement const& uid(
        DiscreteCoordinate<Tags...> const& tuple) noexcept
{
    return tuple.template uid<QueryTag>();
}

template <class QueryTag, class... Tags>
inline constexpr DiscreteCoordinateElement& uid(DiscreteCoordinate<Tags...>& tuple) noexcept
{
    return tuple.template uid<QueryTag>();
}

template <class QueryTag, class... Tags>
inline constexpr DiscreteCoordinateElement const& uid_or(
        DiscreteCoordinate<Tags...> const& tuple,
        DiscreteCoordinateElement const& default_value) noexcept
{
    return tuple.template uid_or<QueryTag>(default_value);
}

template <class... QueryTags, class... Tags>
inline constexpr DiscreteCoordinate<QueryTags...> select(
        DiscreteCoordinate<Tags...> const& arr) noexcept
{
    return DiscreteCoordinate<QueryTags...>(arr);
}

template <class... QueryTags, class... Tags>
inline constexpr DiscreteCoordinate<QueryTags...> select(DiscreteCoordinate<Tags...>&& arr) noexcept
{
    return DiscreteCoordinate<QueryTags...>(std::move(arr));
}

template <class QueryTag, class HeadTag, class... TailTags>
constexpr DiscreteCoordinate<QueryTag> const& take(
        DiscreteCoordinate<HeadTag> const& head,
        DiscreteCoordinate<TailTags> const&... tags)
{
    static_assert(
            !type_seq_contains_v<detail::TypeSeq<HeadTag>, detail::TypeSeq<TailTags...>>,
            "ERROR: tag redundant");
    if constexpr (std::is_same_v<QueryTag, HeadTag>) {
        return head;
    } else {
        static_assert(sizeof...(TailTags) > 0, "ERROR: tag not found");
        return take<QueryTag>(tags...);
    }
}

namespace detail {

/// Returns a reference to the underlying `std::array`
template <class... Tags>
constexpr inline std::array<DiscreteCoordinateElement, sizeof...(Tags)>& array(
        DiscreteCoordinate<Tags...>& v) noexcept
{
    return v.m_values;
}

/// Returns a reference to the underlying `std::array`
template <class... Tags>
constexpr inline std::array<DiscreteCoordinateElement, sizeof...(Tags)> const& array(
        DiscreteCoordinate<Tags...> const& v) noexcept
{
    return v.m_values;
}

} // namespace detail

/** A DiscreteCoordinate identifies an element of the discrete dimension
 *
 * Each one is tagged by its associated dimensions.
 */
template <class... Tags>
class DiscreteCoordinate
{
    using tags_seq = detail::TypeSeq<Tags...>;

    friend constexpr std::array<DiscreteCoordinateElement, sizeof...(Tags)>& detail::array<Tags...>(
            DiscreteCoordinate<Tags...>& v) noexcept;

    friend constexpr std::array<DiscreteCoordinateElement, sizeof...(Tags)> const& detail::array<
            Tags...>(DiscreteCoordinate<Tags...> const& v) noexcept;

private:
    std::array<DiscreteCoordinateElement, sizeof...(Tags)> m_values;

public:
    using value_type = DiscreteCoordinateElement;

    static constexpr std::size_t size() noexcept
    {
        return sizeof...(Tags);
    }

public:
    inline constexpr DiscreteCoordinate() = default;

    inline constexpr DiscreteCoordinate(DiscreteCoordinate const&) = default;

    inline constexpr DiscreteCoordinate(DiscreteCoordinate&&) = default;

    template <class... OTags>
    explicit inline constexpr DiscreteCoordinate(DiscreteCoordinate<OTags> const&... other) noexcept
        : m_values {take<Tags>(other...).uid()...}
    {
    }

    template <class... OTags>
    explicit inline constexpr DiscreteCoordinate(DiscreteCoordinate<OTags...> const& other) noexcept
        : m_values {other.template uid<Tags>()...}
    {
    }

    template <
            class... Params,
            class = std::enable_if_t<(std::is_integral_v<Params> && ...)>,
            class = std::enable_if_t<(!is_discrete_coordinate_v<Params> && ...)>,
            class = std::enable_if_t<sizeof...(Params) == sizeof...(Tags)>>
    explicit inline constexpr DiscreteCoordinate(Params const&... params) noexcept
        : m_values {static_cast<value_type>(params)...}
    {
    }

    constexpr inline DiscreteCoordinate& operator=(DiscreteCoordinate const& other) = default;

    constexpr inline DiscreteCoordinate& operator=(DiscreteCoordinate&& other) = default;

    template <class... OTags>
    constexpr inline DiscreteCoordinate& operator=(
            DiscreteCoordinate<OTags...> const& other) noexcept
    {
        m_values = other.m_values;
        return *this;
    }

    template <class... OTags>
    constexpr inline DiscreteCoordinate& operator=(DiscreteCoordinate<OTags...>&& other) noexcept
    {
        m_values = std::move(other.m_values);
        return *this;
    }

    template <class QueryTag>
    value_type const& uid_or(value_type const& default_value) const&
    {
        if constexpr (in_tags_v<QueryTag, tags_seq>) {
            return m_values[type_seq_rank_v<QueryTag, tags_seq>];
        } else {
            return default_value;
        }
    }

    template <class QueryTag>
    inline constexpr value_type& uid() noexcept
    {
        using namespace detail;
        static_assert(
                in_tags_v<QueryTag, tags_seq>,
                "requested Tag absent from DiscreteCoordinate");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <class QueryTag>
    inline constexpr value_type const& uid() const noexcept
    {
        using namespace detail;
        static_assert(
                in_tags_v<QueryTag, tags_seq>,
                "requested Tag absent from DiscreteCoordinate");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <std::size_t N = sizeof...(Tags)>
    constexpr inline std::enable_if_t<N == 1, value_type&> uid() noexcept
    {
        return m_values[0];
    }

    template <std::size_t N = sizeof...(Tags)>
    constexpr inline std::enable_if_t<N == 1, value_type const&> uid() const noexcept
    {
        return m_values[0];
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    constexpr inline DiscreteCoordinate& operator++()
    {
        ++m_values[0];
        return *this;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    constexpr inline DiscreteCoordinate operator++(int)
    {
        DiscreteCoordinate const tmp = *this;
        ++m_values[0];
        return tmp;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    constexpr inline DiscreteCoordinate& operator--()
    {
        ++m_values[0];
        return *this;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    constexpr inline DiscreteCoordinate operator--(int)
    {
        DiscreteCoordinate const tmp = *this;
        ++m_values[0];
        return tmp;
    }

    template <class... OTags>
    constexpr inline DiscreteCoordinate& operator+=(DiscreteVector<OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, detail::TypeSeq<OTags...>>);
        ((m_values[type_seq_rank_v<Tags, tags_seq>] += rhs.template get<Tags>()), ...);
        return *this;
    }

    template <
            class IntegralType,
            std::size_t N = sizeof...(Tags),
            class = std::enable_if_t<N == 1>,
            class = std::enable_if_t<std::is_integral_v<IntegralType>>>
    constexpr inline DiscreteCoordinate& operator+=(IntegralType const& rhs)
    {
        m_values[0] += rhs;
        return *this;
    }

    template <class... OTags>
    constexpr inline DiscreteCoordinate& operator-=(DiscreteVector<OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, detail::TypeSeq<OTags...>>);
        ((m_values[type_seq_rank_v<Tags, tags_seq>] -= rhs.template get<Tags>()), ...);
        return *this;
    }

    template <
            class IntegralType,
            std::size_t N = sizeof...(Tags),
            class = std::enable_if_t<N == 1>,
            class = std::enable_if_t<std::is_integral_v<IntegralType>>>
    constexpr inline DiscreteCoordinate& operator-=(IntegralType const& rhs)
    {
        m_values[0] -= rhs;
        return *this;
    }
};

inline std::ostream& operator<<(std::ostream& out, DiscreteCoordinate<> const&)
{
    out << "()";
    return out;
}

template <class Head, class... Tags>
std::ostream& operator<<(std::ostream& out, DiscreteCoordinate<Head, Tags...> const& arr)
{
    out << "(";
    out << uid<Head>(arr);
    ((out << ", " << uid<Tags>(arr)), ...);
    out << ")";
    return out;
}


template <class... Tags, class... OTags>
constexpr inline bool operator==(
        DiscreteCoordinate<Tags...> const& lhs,
        DiscreteCoordinate<OTags...> const& rhs) noexcept
{
    return ((lhs.template uid<Tags>() == rhs.template uid<Tags>()) && ...);
}

template <class... Tags, class... OTags>
constexpr inline bool operator!=(
        DiscreteCoordinate<Tags...> const& lhs,
        DiscreteCoordinate<OTags...> const& rhs) noexcept
{
    return !(lhs == rhs);
}

template <class Tag>
constexpr inline bool operator<(
        DiscreteCoordinate<Tag> const& lhs,
        DiscreteCoordinate<Tag> const& rhs)
{
    return lhs.uid() < rhs.uid();
}

template <class Tag>
constexpr inline bool operator<=(
        DiscreteCoordinate<Tag> const& lhs,
        DiscreteCoordinate<Tag> const& rhs)
{
    return lhs.uid() <= rhs.uid();
}

template <class Tag>
constexpr inline bool operator>(
        DiscreteCoordinate<Tag> const& lhs,
        DiscreteCoordinate<Tag> const& rhs)
{
    return lhs.uid() > rhs.uid();
}

template <class Tag>
constexpr inline bool operator>=(
        DiscreteCoordinate<Tag> const& lhs,
        DiscreteCoordinate<Tag> const& rhs)
{
    return lhs.uid() >= rhs.uid();
}

/// right external binary operators: +, -

template <class... Tags, class... OTags>
constexpr inline DiscreteCoordinate<Tags...> operator+(
        DiscreteCoordinate<Tags...> const& lhs,
        DiscreteVector<OTags...> const& rhs)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    return DiscreteCoordinate<Tags...>((uid<Tags>(lhs) + get<Tags>(rhs))...);
}

template <
        class Tag,
        class IntegralType,
        class = std::enable_if_t<std::is_integral_v<IntegralType>>,
        class = std::enable_if_t<!is_discrete_vector_v<IntegralType>>>
constexpr inline DiscreteCoordinate<Tag> operator+(
        DiscreteCoordinate<Tag> const& lhs,
        IntegralType const& rhs)
{
    return DiscreteCoordinate<Tag>(uid<Tag>(lhs) + rhs);
}

template <class... Tags, class... OTags>
constexpr inline DiscreteCoordinate<Tags...> operator-(
        DiscreteCoordinate<Tags...> const& lhs,
        DiscreteVector<OTags...> const& rhs)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    return DiscreteCoordinate<Tags...>((uid<Tags>(lhs) - get<Tags>(rhs))...);
}

template <
        class Tag,
        class IntegralType,
        class = std::enable_if_t<std::is_integral_v<IntegralType>>,
        class = std::enable_if_t<!is_discrete_vector_v<IntegralType>>>
constexpr inline DiscreteCoordinate<Tag> operator-(
        DiscreteCoordinate<Tag> const& lhs,
        IntegralType const& rhs)
{
    return DiscreteCoordinate<Tag>(uid<Tag>(lhs) - rhs);
}

/// binary operator: -

template <class... Tags, class... OTags>
constexpr inline DiscreteVector<Tags...> operator-(
        DiscreteCoordinate<Tags...> const& lhs,
        DiscreteCoordinate<OTags...> const& rhs)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    return DiscreteVector<Tags...>((uid<Tags>(lhs) - uid<Tags>(rhs))...);
}
