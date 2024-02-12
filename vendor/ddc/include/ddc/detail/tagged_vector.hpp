// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstddef>
#include <ostream>
#include <type_traits>
#include <utility>

#include "ddc/detail/macros.hpp"
#include "ddc/detail/type_seq.hpp"

namespace ddc {

namespace detail {

template <class, class...>
class TaggedVector;

template <class T>
struct IsTaggedVector : std::false_type
{
};

template <class ElementType, class... Tags>
struct IsTaggedVector<TaggedVector<ElementType, Tags...>> : std::true_type
{
};

template <class T>
inline constexpr bool is_tagged_vector_v = IsTaggedVector<T>::value;

template <class ElementType, class... Tags>
struct ToTypeSeq<TaggedVector<ElementType, Tags...>>
{
    using type = TypeSeq<Tags...>;
};

} // namespace detail


template <class QueryTag, class ElementType, class... Tags>
KOKKOS_FUNCTION constexpr ElementType const& get(
        detail::TaggedVector<ElementType, Tags...> const& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class ElementType, class... Tags>
KOKKOS_FUNCTION constexpr ElementType& get(
        detail::TaggedVector<ElementType, Tags...>& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class ElementType, class... Tags>
KOKKOS_FUNCTION constexpr ElementType const& get_or(
        detail::TaggedVector<ElementType, Tags...> const& tuple,
        ElementType const& default_value) noexcept
{
    return tuple.template get_or<QueryTag>(default_value);
}

namespace detail {

/// Unary operators: +, -

template <class ElementType, class... Tags>
KOKKOS_FUNCTION constexpr detail::TaggedVector<ElementType, Tags...> operator+(
        detail::TaggedVector<ElementType, Tags...> const& x)
{
    return x;
}

template <class ElementType, class... Tags>
KOKKOS_FUNCTION constexpr detail::TaggedVector<ElementType, Tags...> operator-(
        detail::TaggedVector<ElementType, Tags...> const& x)
{
    return detail::TaggedVector<ElementType, Tags...>((-get<Tags>(x))...);
}

/// Internal binary operators: +, -

template <class ElementType, class... Tags, class OElementType, class... OTags>
KOKKOS_FUNCTION constexpr auto operator+(
        detail::TaggedVector<ElementType, Tags...> const& lhs,
        detail::TaggedVector<OElementType, OTags...> const& rhs)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType = decltype(std::declval<ElementType>() + std::declval<OElementType>());
    return detail::TaggedVector<RElementType, Tags...>((get<Tags>(lhs) + get<Tags>(rhs))...);
}

template <
        class ElementType,
        class Tag,
        class OElementType,
        class = std::enable_if_t<!detail::is_tagged_vector_v<OElementType>>,
        class = std::enable_if_t<std::is_convertible_v<OElementType, ElementType>>>
KOKKOS_FUNCTION constexpr auto operator+(
        detail::TaggedVector<ElementType, Tag> const& lhs,
        OElementType const& rhs)
{
    using RElementType = decltype(std::declval<ElementType>() + std::declval<OElementType>());
    return detail::TaggedVector<RElementType, Tag>(get<Tag>(lhs) + rhs);
}

template <
        class ElementType,
        class Tag,
        class OElementType,
        class = std::enable_if_t<!detail::is_tagged_vector_v<OElementType>>,
        class = std::enable_if_t<std::is_convertible_v<ElementType, OElementType>>>
KOKKOS_FUNCTION constexpr auto operator+(
        OElementType const& lhs,
        detail::TaggedVector<ElementType, Tag> const& rhs)
{
    using RElementType = decltype(std::declval<ElementType>() + std::declval<OElementType>());
    return detail::TaggedVector<RElementType, Tag>(lhs + get<Tag>(rhs));
}

template <class ElementType, class... Tags, class OElementType, class... OTags>
KOKKOS_FUNCTION constexpr auto operator-(
        detail::TaggedVector<ElementType, Tags...> const& lhs,
        detail::TaggedVector<OElementType, OTags...> const& rhs)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType = decltype(std::declval<ElementType>() - std::declval<OElementType>());
    return detail::TaggedVector<RElementType, Tags...>((get<Tags>(lhs) - get<Tags>(rhs))...);
}

template <
        class ElementType,
        class Tag,
        class OElementType,
        class = std::enable_if_t<!detail::is_tagged_vector_v<OElementType>>,
        class = std::enable_if_t<std::is_convertible_v<OElementType, ElementType>>>
KOKKOS_FUNCTION constexpr auto operator-(
        detail::TaggedVector<ElementType, Tag> const& lhs,
        OElementType const& rhs)
{
    using RElementType = decltype(std::declval<ElementType>() + std::declval<OElementType>());
    return detail::TaggedVector<RElementType, Tag>(get<Tag>(lhs) - rhs);
}

template <
        class ElementType,
        class Tag,
        class OElementType,
        class = std::enable_if_t<!detail::is_tagged_vector_v<OElementType>>,
        class = std::enable_if_t<std::is_convertible_v<ElementType, OElementType>>>
KOKKOS_FUNCTION constexpr auto operator-(
        OElementType const& lhs,
        detail::TaggedVector<ElementType, Tag> const& rhs)
{
    using RElementType = decltype(std::declval<ElementType>() + std::declval<OElementType>());
    return detail::TaggedVector<RElementType, Tag>(lhs - get<Tag>(rhs));
}

/// external left binary operator: *

template <
        class ElementType,
        class OElementType,
        class... Tags,
        class = std::enable_if_t<!detail::is_tagged_vector_v<OElementType>>,
        class = std::enable_if_t<std::is_convertible_v<ElementType, OElementType>>>
KOKKOS_FUNCTION constexpr auto operator*(
        ElementType const& lhs,
        detail::TaggedVector<OElementType, Tags...> const& rhs)
{
    using RElementType = decltype(std::declval<ElementType>() * std::declval<OElementType>());
    return detail::TaggedVector<RElementType, Tags...>((lhs * get<Tags>(rhs))...);
}

} // namespace detail

template <class... QueryTags, class ElementType, class... Tags>
KOKKOS_FUNCTION constexpr detail::TaggedVector<ElementType, QueryTags...> select(
        detail::TaggedVector<ElementType, Tags...> const& arr) noexcept
{
    return detail::TaggedVector<ElementType, QueryTags...>(arr);
}

template <class... QueryTags, class ElementType, class... Tags>
KOKKOS_FUNCTION constexpr detail::TaggedVector<ElementType, QueryTags...> select(
        detail::TaggedVector<ElementType, Tags...>&& arr) noexcept
{
    return detail::TaggedVector<ElementType, QueryTags...>(std::move(arr));
}

namespace detail {

/// Returns a reference towards the DiscreteElement that contains the QueryTag
template <
        class QueryTag,
        class HeadTaggedVector,
        class... TailTaggedVectors,
        std::enable_if_t<
                is_tagged_vector_v<
                        HeadTaggedVector> && (is_tagged_vector_v<TailTaggedVectors> && ...),
                int> = 1>
KOKKOS_FUNCTION constexpr auto const& take(
        HeadTaggedVector const& head,
        TailTaggedVectors const&... tail)
{
    DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
    if constexpr (type_seq_contains_v<detail::TypeSeq<QueryTag>, to_type_seq_t<HeadTaggedVector>>) {
        static_assert(
                (!type_seq_contains_v<
                         detail::TypeSeq<QueryTag>,
                         to_type_seq_t<TailTaggedVectors>> && ...),
                "ERROR: tag redundant");
        return head;
    } else {
        static_assert(sizeof...(TailTaggedVectors) > 0, "ERROR: tag not found");
        return take<QueryTag>(tail...);
    }
    DDC_IF_NVCC_THEN_POP
}


template <class T>
class ConversionOperators
{
};

template <class ElementType, class Tag>
class ConversionOperators<TaggedVector<ElementType, Tag>>
{
public:
    KOKKOS_FUNCTION constexpr operator ElementType const &() const noexcept
    {
        return static_cast<TaggedVector<ElementType, Tag> const*>(this)->m_values[0];
    }

    KOKKOS_FUNCTION constexpr operator ElementType&() noexcept
    {
        return static_cast<TaggedVector<ElementType, Tag>*>(this)->m_values[0];
    }
};

template <class ElementType, class... Tags>
class TaggedVector : public ConversionOperators<TaggedVector<ElementType, Tags...>>
{
    friend class ConversionOperators<TaggedVector<ElementType, Tags...>>;

    using tags_seq = detail::TypeSeq<Tags...>;

private:
    std::array<ElementType, sizeof...(Tags)> m_values;

public:
    using value_type = ElementType;

    static KOKKOS_FUNCTION constexpr std::size_t size() noexcept
    {
        return sizeof...(Tags);
    }

public:
    KOKKOS_DEFAULTED_FUNCTION constexpr TaggedVector() = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr TaggedVector(TaggedVector const&) = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr TaggedVector(TaggedVector&&) = default;

    template <class... TVectors, class = std::enable_if_t<(is_tagged_vector_v<TVectors> && ...)>>
    explicit KOKKOS_FUNCTION constexpr TaggedVector(TVectors const&... delems) noexcept
        : m_values {static_cast<ElementType>(take<Tags>(delems...).template get<Tags>())...}
    {
    }

    template <
            class... Params,
            class = std::enable_if_t<(!is_tagged_vector_v<Params> && ...)>,
            class = std::enable_if_t<(std::is_convertible_v<Params, ElementType> && ...)>,
            class = std::enable_if_t<sizeof...(Params) == sizeof...(Tags)>>
    explicit KOKKOS_FUNCTION constexpr TaggedVector(Params const&... params) noexcept
        : m_values {static_cast<ElementType>(params)...}
    {
    }

    KOKKOS_DEFAULTED_FUNCTION ~TaggedVector() = default;

    KOKKOS_DEFAULTED_FUNCTION TaggedVector& operator=(TaggedVector const& other) = default;

    KOKKOS_DEFAULTED_FUNCTION TaggedVector& operator=(TaggedVector&& other) = default;

    template <class... OTags>
    KOKKOS_FUNCTION constexpr TaggedVector& operator=(
            TaggedVector<ElementType, OTags...> const& other) noexcept
    {
        ((this->get<Tags>() = other.template get<Tags>()), ...);
        return *this;
    }

    template <class... OTags>
    KOKKOS_FUNCTION constexpr TaggedVector& operator=(
            TaggedVector<ElementType, OTags...>&& other) noexcept
    {
        ((this->get<Tags>() = std::move(other.template get<Tags>())), ...);
        return *this;
    }

    /// Returns a reference to the underlying `std::array`
    KOKKOS_FUNCTION constexpr std::array<ElementType, sizeof...(Tags)>& array() noexcept
    {
        return m_values;
    }

    /// Returns a const reference to the underlying `std::array`
    KOKKOS_FUNCTION constexpr std::array<ElementType, sizeof...(Tags)> const& array() const noexcept
    {
        return m_values;
    }

    KOKKOS_FUNCTION constexpr ElementType& operator[](size_t pos)
    {
        return m_values[pos];
    }

    KOKKOS_FUNCTION constexpr ElementType const& operator[](size_t pos) const
    {
        return m_values[pos];
    }

    template <class OElementType, class... OTags>
    KOKKOS_FUNCTION constexpr bool operator==(
            TaggedVector<OElementType, OTags...> const& rhs) const noexcept
    {
        return ((m_values[type_seq_rank_v<Tags, tags_seq>] == rhs.template get<Tags>()) && ...);
    }

    template <class OElementType, class... OTags>
    KOKKOS_FUNCTION constexpr bool operator!=(
            TaggedVector<OElementType, OTags...> const& rhs) const noexcept
    {
        return !(*this == rhs);
    }

    template <class QueryTag>
    KOKKOS_FUNCTION constexpr ElementType& get() noexcept
    {
        using namespace detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from TaggedVector");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <class QueryTag>
    KOKKOS_FUNCTION constexpr ElementType const& get() const noexcept
    {
        using namespace detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from TaggedVector");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <class QueryTag>
    KOKKOS_FUNCTION constexpr ElementType const& get_or(ElementType const& default_value) const&
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
    KOKKOS_FUNCTION constexpr std::enable_if_t<N == 1, ElementType const&> value() const noexcept
    {
        return m_values[0];
    }

    template <class OElementType, class... OTags>
    KOKKOS_FUNCTION constexpr TaggedVector& operator+=(
            TaggedVector<OElementType, OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, detail::TypeSeq<OTags...>>);
        ((m_values[type_seq_rank_v<Tags, tags_seq>] += rhs.template get<Tags>()), ...);
        return *this;
    }

    template <class OElementType>
    KOKKOS_FUNCTION constexpr TaggedVector& operator+=(OElementType const& rhs)
    {
        ((m_values[type_seq_rank_v<Tags, tags_seq>] += rhs), ...);
        return *this;
    }

    template <class OElementType, class... OTags>
    KOKKOS_FUNCTION constexpr TaggedVector& operator-=(
            TaggedVector<OElementType, OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, detail::TypeSeq<OTags...>>);
        ((m_values[type_seq_rank_v<Tags, tags_seq>] -= rhs.template get<Tags>()), ...);
        return *this;
    }

    template <class OElementType>
    KOKKOS_FUNCTION constexpr TaggedVector& operator-=(OElementType const& rhs)
    {
        ((m_values[type_seq_rank_v<Tags, tags_seq>] -= rhs), ...);
        return *this;
    }

    template <class OElementType, class... OTags>
    KOKKOS_FUNCTION constexpr TaggedVector& operator*=(
            TaggedVector<OElementType, OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, detail::TypeSeq<OTags...>>);
        ((m_values[type_seq_rank_v<Tags, tags_seq>] *= rhs.template get<Tags>()), ...);
        return *this;
    }
};

template <class ElementType>
std::ostream& operator<<(std::ostream& out, TaggedVector<ElementType> const&)
{
    out << "()";
    return out;
}

template <class ElementType, class Head, class... Tags>
std::ostream& operator<<(std::ostream& out, TaggedVector<ElementType, Head, Tags...> const& arr)
{
    out << "(";
    out << get<Head>(arr);
    ((out << ", " << get<Tags>(arr)), ...);
    out << ")";
    return out;
}

} // namespace detail

} // namespace ddc
