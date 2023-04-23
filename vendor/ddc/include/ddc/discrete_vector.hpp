// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstddef>
#include <ostream>
#include <utility>

#include "ddc/detail/type_seq.hpp"

namespace ddc {

template <class...>
class DiscreteVector;

template <class T>
struct IsDiscreteVector : std::false_type
{
};

template <class... Tags>
struct IsDiscreteVector<DiscreteVector<Tags...>> : std::true_type
{
};

template <class T>
inline constexpr bool is_discrete_vector_v = IsDiscreteVector<T>::value;


/** A DiscreteVectorElement is a scalar that represents the difference between two coordinates.
 */
using DiscreteVectorElement = std::ptrdiff_t;

template <class QueryTag, class... Tags>
inline constexpr DiscreteVectorElement const& get(DiscreteVector<Tags...> const& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class... Tags>
inline constexpr DiscreteVectorElement& get(DiscreteVector<Tags...>& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class... Tags>
inline constexpr DiscreteVectorElement const& get_or(
        DiscreteVector<Tags...> const& tuple,
        DiscreteVectorElement const& default_value) noexcept
{
    return tuple.template get_or<QueryTag>(default_value);
}

/// Unary operators: +, -

template <class... Tags>
constexpr inline DiscreteVector<Tags...> operator+(DiscreteVector<Tags...> const& x)
{
    return x;
}

template <class... Tags>
constexpr inline DiscreteVector<Tags...> operator-(DiscreteVector<Tags...> const& x)
{
    return DiscreteVector<Tags...>((-get<Tags>(x))...);
}

/// Internal binary operators: +, -

template <class... Tags, class... OTags>
constexpr inline DiscreteVector<Tags...> operator+(
        DiscreteVector<Tags...> const& lhs,
        DiscreteVector<OTags...> const& rhs)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    return DiscreteVector<Tags...>((get<Tags>(lhs) + get<Tags>(rhs))...);
}

template <class Tag, class IntegralType, class = std::enable_if_t<std::is_integral_v<IntegralType>>>
constexpr inline DiscreteVector<Tag> operator+(
        DiscreteVector<Tag> const& lhs,
        IntegralType const& rhs)
{
    return DiscreteVector<Tag>(get<Tag>(lhs) + rhs);
}

template <class IntegralType, class Tag, class = std::enable_if_t<std::is_integral_v<IntegralType>>>
constexpr inline DiscreteVector<Tag> operator+(
        IntegralType const& lhs,
        DiscreteVector<Tag> const& rhs)
{
    return DiscreteVector<Tag>(lhs + get<Tag>(rhs));
}

template <class... Tags, class... OTags>
constexpr inline DiscreteVector<Tags...> operator-(
        DiscreteVector<Tags...> const& lhs,
        DiscreteVector<OTags...> const& rhs)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    return DiscreteVector<Tags...>((get<Tags>(lhs) - get<Tags>(rhs))...);
}

template <class Tag, class IntegralType, class = std::enable_if_t<std::is_integral_v<IntegralType>>>
constexpr inline DiscreteVector<Tag> operator-(
        DiscreteVector<Tag> const& lhs,
        IntegralType const& rhs)
{
    return DiscreteVector<Tag>(get<Tag>(lhs) - rhs);
}

template <class IntegralType, class Tag, class = std::enable_if_t<std::is_integral_v<IntegralType>>>
constexpr inline DiscreteVector<Tag> operator-(
        IntegralType const& lhs,
        DiscreteVector<Tag> const& rhs)
{
    return DiscreteVector<Tag>(lhs - get<Tag>(rhs));
}

/// external left binary operator: *

template <
        class IntegralType,
        class... Tags,
        class = std::enable_if_t<std::is_integral_v<IntegralType>>>
constexpr inline auto operator*(IntegralType const& lhs, DiscreteVector<Tags...> const& rhs)
{
    return DiscreteVector<Tags...>((lhs * get<Tags>(rhs))...);
}

template <class... QueryTags, class... Tags>
inline constexpr DiscreteVector<QueryTags...> select(DiscreteVector<Tags...> const& arr) noexcept
{
    return DiscreteVector<QueryTags...>(arr);
}

template <class... QueryTags, class... Tags>
inline constexpr DiscreteVector<QueryTags...> select(DiscreteVector<Tags...>&& arr) noexcept
{
    return DiscreteVector<QueryTags...>(std::move(arr));
}

template <class QueryTag, class HeadTag, class... TailTags>
constexpr DiscreteVector<QueryTag> const& take(
        DiscreteVector<HeadTag> const& head,
        DiscreteVector<TailTags> const&... tags)
{
    DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
    static_assert(
            !type_seq_contains_v<detail::TypeSeq<HeadTag>, detail::TypeSeq<TailTags...>>,
            "ERROR: tag redundant");
    if constexpr (std::is_same_v<QueryTag, HeadTag>) {
        return head;
    } else {
        static_assert(sizeof...(TailTags) > 0, "ERROR: tag not found");
        return take<QueryTag>(tags...);
    }
    DDC_IF_NVCC_THEN_POP
}

template <class T>
class ConversionOperators
{
};

template <class Tag>
class ConversionOperators<DiscreteVector<Tag>>
{
public:
    constexpr inline operator DiscreteVectorElement const &() const noexcept
    {
        return static_cast<DiscreteVector<Tag> const*>(this)->m_values[0];
    }

    constexpr inline operator DiscreteVectorElement&() noexcept
    {
        return static_cast<DiscreteVector<Tag>*>(this)->m_values[0];
    }
};

namespace detail {

/// Returns a reference to the underlying `std::array`
template <class... Tags>
constexpr inline std::array<DiscreteVectorElement, sizeof...(Tags)>& array(
        DiscreteVector<Tags...>& v) noexcept
{
    return v.m_values;
}

/// Returns a reference to the underlying `std::array`
template <class... Tags>
constexpr inline std::array<DiscreteVectorElement, sizeof...(Tags)> const& array(
        DiscreteVector<Tags...> const& v) noexcept
{
    return v.m_values;
}

} // namespace detail

/** A DiscreteVector is a vector in the discrete dimension
 *
 * Each is tagged by its associated dimensions.
 */
template <class... Tags>
class DiscreteVector : public ConversionOperators<DiscreteVector<Tags...>>
{
    friend class ConversionOperators<DiscreteVector<Tags...>>;

    friend constexpr std::array<DiscreteVectorElement, sizeof...(Tags)>& detail::array<Tags...>(
            DiscreteVector<Tags...>& v) noexcept;
    friend constexpr std::array<DiscreteVectorElement, sizeof...(Tags)> const& detail::array<
            Tags...>(DiscreteVector<Tags...> const& v) noexcept;

    using tags_seq = detail::TypeSeq<Tags...>;

private:
    std::array<DiscreteVectorElement, sizeof...(Tags)> m_values;

public:
    static constexpr std::size_t size() noexcept
    {
        return sizeof...(Tags);
    }

public:
    inline constexpr DiscreteVector() = default;

    inline constexpr DiscreteVector(DiscreteVector const&) = default;

    inline constexpr DiscreteVector(DiscreteVector&&) = default;

    template <class... OTags>
    explicit inline constexpr DiscreteVector(DiscreteVector<OTags> const&... other) noexcept
        : m_values {take<Tags>(other...).value()...}
    {
    }

    template <class... OTags>
    explicit inline constexpr DiscreteVector(DiscreteVector<OTags...> const& other) noexcept
        : m_values {other.template get<Tags>()...}
    {
    }

    template <
            class... Params,
            class = std::enable_if_t<(std::is_convertible_v<Params, DiscreteVectorElement> && ...)>,
            class = std::enable_if_t<(!is_discrete_vector_v<Params> && ...)>,
            class = std::enable_if_t<sizeof...(Params) == sizeof...(Tags)>>
    explicit inline constexpr DiscreteVector(Params const&... params) noexcept
        : m_values {static_cast<DiscreteVectorElement>(params)...}
    {
    }

    constexpr inline DiscreteVector& operator=(DiscreteVector const& other) = default;

    constexpr inline DiscreteVector& operator=(DiscreteVector&& other) = default;

    template <class... OTags>
    constexpr inline DiscreteVector& operator=(DiscreteVector<OTags...> const& other) noexcept
    {
        m_values = other.m_values;
        return *this;
    }

    template <class... OTags>
    constexpr inline DiscreteVector& operator=(DiscreteVector<OTags...>&& other) noexcept
    {
        m_values = std::move(other.m_values);
        return *this;
    }

    template <class... OTags>
    constexpr inline bool operator==(DiscreteVector<OTags...> const& rhs) const noexcept
    {
        return ((m_values[type_seq_rank_v<Tags, tags_seq>] == rhs.template get<Tags>()) && ...);
    }

    template <class... OTags>
    constexpr inline bool operator!=(DiscreteVector<OTags...> const& rhs) const noexcept
    {
        return !(*this == rhs);
    }

    template <class QueryTag>
    inline constexpr DiscreteVectorElement& get() noexcept
    {
        using namespace detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from DiscreteVector");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <class QueryTag>
    inline constexpr DiscreteVectorElement const& get() const noexcept
    {
        using namespace detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from DiscreteVector");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <class QueryTag>
    DiscreteVectorElement const& get_or(DiscreteVectorElement const& default_value) const&
    {
        DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
        if constexpr (in_tags_v<QueryTag, tags_seq>) {
            return m_values[type_seq_rank_v<QueryTag, tags_seq>];
        } else {
            return default_value;
        }
        DDC_IF_NVCC_THEN_POP
    }

    template <std::size_t N = sizeof...(Tags)>
    constexpr inline std::enable_if_t<N == 1, DiscreteVectorElement const&> value() const noexcept
    {
        return m_values[0];
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    constexpr inline DiscreteVector& operator++()
    {
        ++m_values[0];
        return *this;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    constexpr inline DiscreteVector operator++(int)
    {
        DiscreteVector const tmp = *this;
        ++m_values[0];
        return tmp;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    constexpr inline DiscreteVector& operator--()
    {
        ++m_values[0];
        return *this;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    constexpr inline DiscreteVector operator--(int)
    {
        DiscreteVector const tmp = *this;
        ++m_values[0];
        return tmp;
    }

    template <class... OTags>
    constexpr inline DiscreteVector& operator+=(DiscreteVector<OTags...> const& rhs)
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
    constexpr inline DiscreteVector& operator+=(IntegralType const& rhs)
    {
        m_values[0] += rhs;
        return *this;
    }

    template <class... OTags>
    constexpr inline DiscreteVector& operator-=(DiscreteVector<OTags...> const& rhs)
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
    constexpr inline DiscreteVector& operator-=(IntegralType const& rhs)
    {
        m_values[0] -= rhs;
        return *this;
    }

    template <class... OTags>
    constexpr inline DiscreteVector& operator*=(DiscreteVector<OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, detail::TypeSeq<OTags...>>);
        ((m_values[type_seq_rank_v<Tags, tags_seq>] *= rhs.template get<Tags>()), ...);
        return *this;
    }
};

template <class Tag>
constexpr inline bool operator<(DiscreteVector<Tag> const& lhs, DiscreteVector<Tag> const& rhs)
{
    return lhs.value() < rhs.value();
}

template <class Tag, class IntegralType>
constexpr inline bool operator<(DiscreteVector<Tag> const& lhs, IntegralType const& rhs)
{
    return lhs.value() < rhs;
}

inline std::ostream& operator<<(std::ostream& out, DiscreteVector<> const&)
{
    out << "()";
    return out;
}

template <class Head, class... Tags>
std::ostream& operator<<(std::ostream& out, DiscreteVector<Head, Tags...> const& arr)
{
    out << "(";
    out << get<Head>(arr);
    ((out << ", " << get<Tags>(arr)), ...);
    out << ")";
    return out;
}
} // namespace ddc
