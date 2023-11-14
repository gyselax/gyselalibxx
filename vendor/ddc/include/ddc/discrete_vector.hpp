// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstddef>
#include <ostream>
#include <utility>

#include <Kokkos_Core.hpp>

#include "ddc/detail/macros.hpp"
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


namespace detail {

template <class... Tags>
struct ToTypeSeq<DiscreteVector<Tags...>>
{
    using type = TypeSeq<Tags...>;
};

} // namespace detail

/** A DiscreteVectorElement is a scalar that represents the difference between two coordinates.
 */
using DiscreteVectorElement = std::ptrdiff_t;

template <class QueryTag, class... Tags>
KOKKOS_FUNCTION constexpr DiscreteVectorElement const& get(
        DiscreteVector<Tags...> const& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class... Tags>
KOKKOS_FUNCTION constexpr DiscreteVectorElement& get(DiscreteVector<Tags...>& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class... Tags>
KOKKOS_FUNCTION constexpr DiscreteVectorElement const& get_or(
        DiscreteVector<Tags...> const& tuple,
        DiscreteVectorElement const& default_value) noexcept
{
    return tuple.template get_or<QueryTag>(default_value);
}

/// Unary operators: +, -

template <class... Tags>
KOKKOS_FUNCTION constexpr DiscreteVector<Tags...> operator+(DiscreteVector<Tags...> const& x)
{
    return x;
}

template <class... Tags>
KOKKOS_FUNCTION constexpr DiscreteVector<Tags...> operator-(DiscreteVector<Tags...> const& x)
{
    return DiscreteVector<Tags...>((-get<Tags>(x))...);
}

/// Internal binary operators: +, -

template <class... Tags, class... OTags>
KOKKOS_FUNCTION constexpr auto operator+(
        DiscreteVector<Tags...> const& lhs,
        DiscreteVector<OTags...> const& rhs)
{
    using detail::TypeSeq;
    if constexpr (sizeof...(Tags) >= sizeof...(OTags)) {
        static_assert(((type_seq_contains_v<TypeSeq<OTags>, TypeSeq<Tags...>>)&&...));
        DiscreteVector<Tags...> result(lhs);
        result += rhs;
        return result;
    } else {
        static_assert(((type_seq_contains_v<TypeSeq<Tags>, TypeSeq<OTags...>>)&&...));
        DiscreteVector<OTags...> result(rhs);
        result += lhs;
        return result;
    }
}

template <class... Tags, class... OTags>
KOKKOS_FUNCTION constexpr auto operator-(
        DiscreteVector<Tags...> const& lhs,
        DiscreteVector<OTags...> const& rhs)
{
    using detail::TypeSeq;
    if constexpr (sizeof...(Tags) >= sizeof...(OTags)) {
        static_assert(((type_seq_contains_v<TypeSeq<OTags>, TypeSeq<Tags...>>)&&...));
        DiscreteVector<Tags...> result(lhs);
        result -= rhs;
        return result;
    } else {
        static_assert(((type_seq_contains_v<TypeSeq<Tags>, TypeSeq<OTags...>>)&&...));
        DiscreteVector<OTags...> result(-rhs);
        result += lhs;
        return result;
    }
}

template <class Tag, class IntegralType, class = std::enable_if_t<std::is_integral_v<IntegralType>>>
KOKKOS_FUNCTION constexpr DiscreteVector<Tag> operator+(
        DiscreteVector<Tag> const& lhs,
        IntegralType const& rhs)
{
    return DiscreteVector<Tag>(get<Tag>(lhs) + rhs);
}

template <class IntegralType, class Tag, class = std::enable_if_t<std::is_integral_v<IntegralType>>>
KOKKOS_FUNCTION constexpr DiscreteVector<Tag> operator+(
        IntegralType const& lhs,
        DiscreteVector<Tag> const& rhs)
{
    return DiscreteVector<Tag>(lhs + get<Tag>(rhs));
}

template <class Tag, class IntegralType, class = std::enable_if_t<std::is_integral_v<IntegralType>>>
KOKKOS_FUNCTION constexpr DiscreteVector<Tag> operator-(
        DiscreteVector<Tag> const& lhs,
        IntegralType const& rhs)
{
    return DiscreteVector<Tag>(get<Tag>(lhs) - rhs);
}

template <class IntegralType, class Tag, class = std::enable_if_t<std::is_integral_v<IntegralType>>>
KOKKOS_FUNCTION constexpr DiscreteVector<Tag> operator-(
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
KOKKOS_FUNCTION constexpr auto operator*(
        IntegralType const& lhs,
        DiscreteVector<Tags...> const& rhs)
{
    return DiscreteVector<Tags...>((lhs * get<Tags>(rhs))...);
}

template <class... QueryTags, class... Tags>
KOKKOS_FUNCTION constexpr DiscreteVector<QueryTags...> select(
        DiscreteVector<Tags...> const& arr) noexcept
{
    return DiscreteVector<QueryTags...>(arr);
}

template <class... QueryTags, class... Tags>
KOKKOS_FUNCTION constexpr DiscreteVector<QueryTags...> select(
        DiscreteVector<Tags...>&& arr) noexcept
{
    return DiscreteVector<QueryTags...>(std::move(arr));
}

/// Returns a reference towards the DiscreteVector that contains the QueryTag
template <
        class QueryTag,
        class HeadDVect,
        class... TailDVects,
        std::enable_if_t<
                is_discrete_vector_v<HeadDVect> && (is_discrete_vector_v<TailDVects> && ...),
                int> = 1>
KOKKOS_FUNCTION constexpr auto const& take(HeadDVect const& head, TailDVects const&... tail)
{
    DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
    if constexpr (type_seq_contains_v<detail::TypeSeq<QueryTag>, to_type_seq_t<HeadDVect>>) {
        static_assert(
                (!type_seq_contains_v<detail::TypeSeq<QueryTag>, to_type_seq_t<TailDVects>> && ...),
                "ERROR: tag redundant");
        return head;
    } else {
        static_assert(sizeof...(TailDVects) > 0, "ERROR: tag not found");
        return take<QueryTag>(tail...);
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
    KOKKOS_FUNCTION constexpr operator DiscreteVectorElement const &() const noexcept
    {
        return static_cast<DiscreteVector<Tag> const*>(this)->m_values[0];
    }

    KOKKOS_FUNCTION constexpr operator DiscreteVectorElement&() noexcept
    {
        return static_cast<DiscreteVector<Tag>*>(this)->m_values[0];
    }
};

namespace detail {

/// Returns a reference to the underlying `std::array`
template <class... Tags>
KOKKOS_FUNCTION constexpr std::array<DiscreteVectorElement, sizeof...(Tags)>& array(
        DiscreteVector<Tags...>& v) noexcept
{
    return v.m_values;
}

/// Returns a reference to the underlying `std::array`
template <class... Tags>
KOKKOS_FUNCTION constexpr std::array<DiscreteVectorElement, sizeof...(Tags)> const& array(
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

    friend KOKKOS_FUNCTION constexpr std::array<DiscreteVectorElement, sizeof...(Tags)>& detail::
            array<Tags...>(DiscreteVector<Tags...>& v) noexcept;
    friend KOKKOS_FUNCTION constexpr std::array<DiscreteVectorElement, sizeof...(Tags)> const&
    detail::array<Tags...>(DiscreteVector<Tags...> const& v) noexcept;

    using tags_seq = detail::TypeSeq<Tags...>;

private:
    std::array<DiscreteVectorElement, sizeof...(Tags)> m_values;

public:
    static KOKKOS_FUNCTION constexpr std::size_t size() noexcept
    {
        return sizeof...(Tags);
    }

public:
    KOKKOS_DEFAULTED_FUNCTION constexpr DiscreteVector() = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr DiscreteVector(DiscreteVector const&) = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr DiscreteVector(DiscreteVector&&) = default;

    template <class... DVects, class = std::enable_if_t<(is_discrete_vector_v<DVects> && ...)>>
    explicit KOKKOS_FUNCTION constexpr DiscreteVector(DVects const&... delems) noexcept
        : m_values {take<Tags>(delems...).template get<Tags>()...}
    {
    }

    template <
            class... Params,
            class = std::enable_if_t<(!is_discrete_vector_v<Params> && ...)>,
            class = std::enable_if_t<(std::is_convertible_v<Params, DiscreteVectorElement> && ...)>,
            class = std::enable_if_t<sizeof...(Params) == sizeof...(Tags)>>
    explicit KOKKOS_FUNCTION constexpr DiscreteVector(Params const&... params) noexcept
        : m_values {static_cast<DiscreteVectorElement>(params)...}
    {
    }

    KOKKOS_DEFAULTED_FUNCTION ~DiscreteVector() = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr DiscreteVector& operator=(DiscreteVector const& other)
            = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr DiscreteVector& operator=(DiscreteVector&& other) = default;

    template <class... OTags>
    KOKKOS_FUNCTION constexpr DiscreteVector& operator=(
            DiscreteVector<OTags...> const& other) noexcept
    {
        m_values = other.m_values;
        return *this;
    }

    template <class... OTags>
    KOKKOS_FUNCTION constexpr DiscreteVector& operator=(DiscreteVector<OTags...>&& other) noexcept
    {
        m_values = std::move(other.m_values);
        return *this;
    }

    template <class... OTags>
    KOKKOS_FUNCTION constexpr bool operator==(DiscreteVector<OTags...> const& rhs) const noexcept
    {
        return ((m_values[type_seq_rank_v<Tags, tags_seq>] == rhs.template get<Tags>()) && ...);
    }

    template <class... OTags>
    KOKKOS_FUNCTION constexpr bool operator!=(DiscreteVector<OTags...> const& rhs) const noexcept
    {
        return !(*this == rhs);
    }

    template <class QueryTag>
    KOKKOS_FUNCTION constexpr DiscreteVectorElement& get() noexcept
    {
        using namespace detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from DiscreteVector");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <class QueryTag>
    KOKKOS_FUNCTION constexpr DiscreteVectorElement const& get() const noexcept
    {
        using namespace detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from DiscreteVector");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <class QueryTag>
    KOKKOS_FUNCTION constexpr DiscreteVectorElement const& get_or(
            DiscreteVectorElement const& default_value) const&
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
    KOKKOS_FUNCTION constexpr std::enable_if_t<N == 1, DiscreteVectorElement const&> value()
            const noexcept
    {
        return m_values[0];
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    KOKKOS_FUNCTION constexpr DiscreteVector& operator++()
    {
        ++m_values[0];
        return *this;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    KOKKOS_FUNCTION constexpr DiscreteVector operator++(int)
    {
        DiscreteVector const tmp = *this;
        ++m_values[0];
        return tmp;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    KOKKOS_FUNCTION constexpr DiscreteVector& operator--()
    {
        ++m_values[0];
        return *this;
    }

    template <std::size_t N = sizeof...(Tags), class = std::enable_if_t<N == 1>>
    KOKKOS_FUNCTION constexpr DiscreteVector operator--(int)
    {
        DiscreteVector const tmp = *this;
        ++m_values[0];
        return tmp;
    }

    template <class... OTags>
    KOKKOS_FUNCTION constexpr DiscreteVector& operator+=(DiscreteVector<OTags...> const& rhs)
    {
        static_assert(((type_seq_contains_v<detail::TypeSeq<OTags>, tags_seq>)&&...));
        ((m_values[type_seq_rank_v<OTags, tags_seq>] += rhs.template get<OTags>()), ...);
        return *this;
    }

    template <
            class IntegralType,
            std::size_t N = sizeof...(Tags),
            class = std::enable_if_t<N == 1>,
            class = std::enable_if_t<std::is_integral_v<IntegralType>>>
    KOKKOS_FUNCTION constexpr DiscreteVector& operator+=(IntegralType const& rhs)
    {
        m_values[0] += rhs;
        return *this;
    }

    template <class... OTags>
    KOKKOS_FUNCTION constexpr DiscreteVector& operator-=(DiscreteVector<OTags...> const& rhs)
    {
        static_assert(((type_seq_contains_v<detail::TypeSeq<OTags>, tags_seq>)&&...));
        ((m_values[type_seq_rank_v<OTags, tags_seq>] -= rhs.template get<OTags>()), ...);
        return *this;
    }

    template <
            class IntegralType,
            std::size_t N = sizeof...(Tags),
            class = std::enable_if_t<N == 1>,
            class = std::enable_if_t<std::is_integral_v<IntegralType>>>
    KOKKOS_FUNCTION constexpr DiscreteVector& operator-=(IntegralType const& rhs)
    {
        m_values[0] -= rhs;
        return *this;
    }

    template <class... OTags>
    KOKKOS_FUNCTION constexpr DiscreteVector& operator*=(DiscreteVector<OTags...> const& rhs)
    {
        static_assert(((type_seq_contains_v<detail::TypeSeq<OTags>, tags_seq>)&&...));
        ((m_values[type_seq_rank_v<OTags, tags_seq>] *= rhs.template get<OTags>()), ...);
        return *this;
    }
};

template <class Tag>
KOKKOS_FUNCTION constexpr bool operator<(
        DiscreteVector<Tag> const& lhs,
        DiscreteVector<Tag> const& rhs)
{
    return lhs.value() < rhs.value();
}

template <class Tag, class IntegralType>
KOKKOS_FUNCTION constexpr bool operator<(DiscreteVector<Tag> const& lhs, IntegralType const& rhs)
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
