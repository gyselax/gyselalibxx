// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstddef>
#include <ostream>
#include <utility>

#include "ddc/detail/macros.hpp"
#include "ddc/detail/type_seq.hpp"

namespace ddc {

namespace ddc_detail {

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


} // namespace ddc_detail


template <class QueryTag, class ElementType, class... Tags>
inline constexpr ElementType const& get(
        ddc_detail::TaggedVector<ElementType, Tags...> const& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class ElementType, class... Tags>
inline constexpr ElementType& get(ddc_detail::TaggedVector<ElementType, Tags...>& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class ElementType, class... Tags>
inline constexpr ElementType const& get_or(
        ddc_detail::TaggedVector<ElementType, Tags...> const& tuple,
        ElementType const& default_value) noexcept
{
    return tuple.template get_or<QueryTag>(default_value);
}

/// Unary operators: +, -

template <class ElementType, class... Tags>
constexpr inline ddc_detail::TaggedVector<ElementType, Tags...> operator+(
        ddc_detail::TaggedVector<ElementType, Tags...> const& x)
{
    return x;
}

template <class ElementType, class... Tags>
constexpr inline ddc_detail::TaggedVector<ElementType, Tags...> operator-(
        ddc_detail::TaggedVector<ElementType, Tags...> const& x)
{
    return ddc_detail::TaggedVector<ElementType, Tags...>((-get<Tags>(x))...);
}

/// Internal binary operators: +, -

template <class ElementType, class... Tags, class OElementType, class... OTags>
constexpr inline auto operator+(
        ddc_detail::TaggedVector<ElementType, Tags...> const& lhs,
        ddc_detail::TaggedVector<OElementType, OTags...> const& rhs)
{
    static_assert(type_seq_same_v<ddc_detail::TypeSeq<Tags...>, ddc_detail::TypeSeq<OTags...>>);
    using RElementType = decltype(std::declval<ElementType>() + std::declval<OElementType>());
    return ddc_detail::TaggedVector<RElementType, Tags...>((get<Tags>(lhs) + get<Tags>(rhs))...);
}

template <
        class ElementType,
        class Tag,
        class OElementType,
        class = std::enable_if_t<!ddc_detail::is_tagged_vector_v<OElementType>>,
        class = std::enable_if_t<std::is_convertible_v<OElementType, ElementType>>>
constexpr inline auto operator+(
        ddc_detail::TaggedVector<ElementType, Tag> const& lhs,
        OElementType const& rhs)
{
    using RElementType = decltype(std::declval<ElementType>() + std::declval<OElementType>());
    return ddc_detail::TaggedVector<RElementType, Tag>(get<Tag>(lhs) + rhs);
}

template <
        class ElementType,
        class Tag,
        class OElementType,
        class = std::enable_if_t<!ddc_detail::is_tagged_vector_v<OElementType>>,
        class = std::enable_if_t<std::is_convertible_v<ElementType, OElementType>>>
constexpr inline auto operator+(
        OElementType const& lhs,
        ddc_detail::TaggedVector<ElementType, Tag> const& rhs)
{
    using RElementType = decltype(std::declval<ElementType>() + std::declval<OElementType>());
    return ddc_detail::TaggedVector<RElementType, Tag>(lhs + get<Tag>(rhs));
}

template <class ElementType, class... Tags, class OElementType, class... OTags>
constexpr inline auto operator-(
        ddc_detail::TaggedVector<ElementType, Tags...> const& lhs,
        ddc_detail::TaggedVector<OElementType, OTags...> const& rhs)
{
    static_assert(type_seq_same_v<ddc_detail::TypeSeq<Tags...>, ddc_detail::TypeSeq<OTags...>>);
    using RElementType = decltype(std::declval<ElementType>() - std::declval<OElementType>());
    return ddc_detail::TaggedVector<RElementType, Tags...>((get<Tags>(lhs) - get<Tags>(rhs))...);
}

template <
        class ElementType,
        class Tag,
        class OElementType,
        class = std::enable_if_t<!ddc_detail::is_tagged_vector_v<OElementType>>,
        class = std::enable_if_t<std::is_convertible_v<OElementType, ElementType>>>
constexpr inline auto operator-(
        ddc_detail::TaggedVector<ElementType, Tag> const& lhs,
        OElementType const& rhs)
{
    using RElementType = decltype(std::declval<ElementType>() + std::declval<OElementType>());
    return ddc_detail::TaggedVector<RElementType, Tag>(get<Tag>(lhs) - rhs);
}

template <
        class ElementType,
        class Tag,
        class OElementType,
        class = std::enable_if_t<!ddc_detail::is_tagged_vector_v<OElementType>>,
        class = std::enable_if_t<std::is_convertible_v<ElementType, OElementType>>>
constexpr inline auto operator-(
        OElementType const& lhs,
        ddc_detail::TaggedVector<ElementType, Tag> const& rhs)
{
    using RElementType = decltype(std::declval<ElementType>() + std::declval<OElementType>());
    return ddc_detail::TaggedVector<RElementType, Tag>(lhs - get<Tag>(rhs));
}

/// external left binary operator: *

template <
        class ElementType,
        class OElementType,
        class... Tags,
        class = std::enable_if_t<!ddc_detail::is_tagged_vector_v<OElementType>>,
        class = std::enable_if_t<std::is_convertible_v<ElementType, OElementType>>>
constexpr inline auto operator*(
        ElementType const& lhs,
        ddc_detail::TaggedVector<OElementType, Tags...> const& rhs)
{
    using RElementType = decltype(std::declval<ElementType>() * std::declval<OElementType>());
    return ddc_detail::TaggedVector<RElementType, Tags...>((lhs * get<Tags>(rhs))...);
}

template <class... QueryTags, class ElementType, class... Tags>
inline constexpr ddc_detail::TaggedVector<ElementType, QueryTags...> select(
        ddc_detail::TaggedVector<ElementType, Tags...> const& arr) noexcept
{
    return ddc_detail::TaggedVector<ElementType, QueryTags...>(arr);
}

template <class... QueryTags, class ElementType, class... Tags>
inline constexpr ddc_detail::TaggedVector<ElementType, QueryTags...> select(
        ddc_detail::TaggedVector<ElementType, Tags...>&& arr) noexcept
{
    return ddc_detail::TaggedVector<ElementType, QueryTags...>(std::move(arr));
}

template <class QueryTag, class ElementType, class HeadTag, class... TailTags>
constexpr ddc_detail::TaggedVector<ElementType, QueryTag> const& take_impl(
        ddc_detail::TaggedVector<ElementType, HeadTag> const& head,
        ddc_detail::TaggedVector<ElementType, TailTags> const&... tags)
{
    DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
    static_assert(
            !type_seq_contains_v<ddc_detail::TypeSeq<HeadTag>, ddc_detail::TypeSeq<TailTags...>>,
            "ERROR: tag redundant");
    if constexpr (std::is_same_v<QueryTag, HeadTag>) {
        return head;
    } else {
        static_assert(sizeof...(TailTags) > 0, "ERROR: tag not found");
        return take_impl<QueryTag>(tags...);
    }
    DDC_IF_NVCC_THEN_POP
}

template <class QueryTag, class ElementType, class... Tags>
constexpr ddc_detail::TaggedVector<ElementType, QueryTag> const& take(
        ddc_detail::TaggedVector<ElementType, Tags> const&... tags)
{
    return take_impl<QueryTag>(tags...);
}


namespace ddc_detail {


template <class T>
class ConversionOperators
{
};

template <class ElementType, class Tag>
class ConversionOperators<TaggedVector<ElementType, Tag>>
{
public:
    constexpr inline operator ElementType const &() const noexcept
    {
        return static_cast<TaggedVector<ElementType, Tag> const*>(this)->m_values[0];
    }

    constexpr inline operator ElementType&() noexcept
    {
        return static_cast<TaggedVector<ElementType, Tag>*>(this)->m_values[0];
    }
};

template <class ElementType, class... Tags>
class TaggedVector : public ConversionOperators<TaggedVector<ElementType, Tags...>>
{
    friend class ConversionOperators<TaggedVector<ElementType, Tags...>>;

    using tags_seq = ddc_detail::TypeSeq<Tags...>;

private:
    std::array<ElementType, sizeof...(Tags)> m_values;

public:
    static constexpr std::size_t size() noexcept
    {
        return sizeof...(Tags);
    }

public:
    inline constexpr TaggedVector() = default;

    inline constexpr TaggedVector(TaggedVector const&) = default;

    inline constexpr TaggedVector(TaggedVector&&) = default;

    template <class... OTags>
    inline constexpr TaggedVector(TaggedVector<ElementType, OTags> const&... other) noexcept
        : m_values {take<Tags>(other...).value()...}
    {
    }

    template <class OElementType, class... OTags>
    explicit inline constexpr TaggedVector(
            TaggedVector<OElementType, OTags...> const& other) noexcept
        : m_values {(static_cast<ElementType>(other.template get<Tags>()))...}
    {
    }

    template <
            class... Params,
            class = std::enable_if_t<(std::is_convertible_v<Params, ElementType> && ...)>,
            class = std::enable_if_t<(!is_tagged_vector_v<Params> && ...)>,
            class = std::enable_if_t<sizeof...(Params) == sizeof...(Tags)>>
    explicit inline constexpr TaggedVector(Params const&... params) noexcept
        : m_values {static_cast<ElementType>(params)...}
    {
    }

    constexpr inline TaggedVector& operator=(TaggedVector const& other) = default;

    constexpr inline TaggedVector& operator=(TaggedVector&& other) = default;

    template <class... OTags>
    constexpr inline TaggedVector& operator=(
            TaggedVector<ElementType, OTags...> const& other) noexcept
    {
        m_values = other.m_values;
        return *this;
    }

    template <class... OTags>
    constexpr inline TaggedVector& operator=(TaggedVector<ElementType, OTags...>&& other) noexcept
    {
        m_values = std::move(other.m_values);
        return *this;
    }

    /// Returns a reference to the underlying `std::array`
    constexpr inline std::array<ElementType, sizeof...(Tags)>& array() noexcept
    {
        return m_values;
    }

    /// Returns a const reference to the underlying `std::array`
    constexpr inline std::array<ElementType, sizeof...(Tags)> const& array() const noexcept
    {
        return m_values;
    }

    constexpr inline ElementType& operator[](size_t pos)
    {
        return m_values[pos];
    }

    constexpr inline ElementType const& operator[](size_t pos) const
    {
        return m_values[pos];
    }

    template <class OElementType, class... OTags>
    constexpr inline bool operator==(TaggedVector<OElementType, OTags...> const& rhs) const noexcept
    {
        return ((m_values[type_seq_rank_v<Tags, tags_seq>] == rhs.template get<Tags>()) && ...);
    }

    template <class OElementType, class... OTags>
    constexpr inline bool operator!=(TaggedVector<OElementType, OTags...> const& rhs) const noexcept
    {
        return !(*this == rhs);
    }

    template <class QueryTag>
    inline constexpr ElementType& get() noexcept
    {
        using namespace ddc_detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from TaggedVector");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <class QueryTag>
    inline constexpr ElementType const& get() const noexcept
    {
        using namespace ddc_detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from TaggedVector");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <class QueryTag>
    ElementType const& get_or(ElementType const& default_value) const&
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
    constexpr inline std::enable_if_t<N == 1, ElementType const&> value() const noexcept
    {
        return m_values[0];
    }

    template <class OElementType, class... OTags>
    constexpr inline TaggedVector& operator+=(TaggedVector<OElementType, OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, ddc_detail::TypeSeq<OTags...>>);
        ((m_values[type_seq_rank_v<Tags, tags_seq>] += rhs.template get<Tags>()), ...);
        return *this;
    }

    template <class OElementType>
    constexpr inline TaggedVector& operator+=(OElementType const& rhs)
    {
        ((m_values[type_seq_rank_v<Tags, tags_seq>] += rhs), ...);
        return *this;
    }

    template <class OElementType, class... OTags>
    constexpr inline TaggedVector& operator-=(TaggedVector<OElementType, OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, ddc_detail::TypeSeq<OTags...>>);
        ((m_values[type_seq_rank_v<Tags, tags_seq>] -= rhs.template get<Tags>()), ...);
        return *this;
    }

    template <class OElementType>
    constexpr inline TaggedVector& operator-=(OElementType const& rhs)
    {
        ((m_values[type_seq_rank_v<Tags, tags_seq>] -= rhs), ...);
        return *this;
    }

    template <class OElementType, class... OTags>
    constexpr inline TaggedVector& operator*=(TaggedVector<OElementType, OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, ddc_detail::TypeSeq<OTags...>>);
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

} // namespace ddc_detail

} // namespace ddc
