#pragma once

#include <array>
#include <cstddef>
#include <ostream>
#include <utility>

#include "ddc/type_seq.hpp"

template <class, class...>
class TaggedVector;

template <class QueryTag, class ElementType, class HeadTag, class... TailTags>
TaggedVector<ElementType, QueryTag> const& take_first(
        TaggedVector<ElementType, HeadTag> const& head,
        TaggedVector<ElementType, TailTags> const&... tags)
{
    if constexpr (std::is_same_v<QueryTag, HeadTag>) {
        return head;
    } else {
        return take_first<QueryTag>(tags...);
    }
}

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

    using tags_seq = detail::TypeSeq<Tags...>;

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

    // template <class... OTags>
    // inline constexpr TaggedVector(TaggedVector<ElementType, OTags> const&... other) noexcept
    //     : m_values {take_first<Tags>(other).value()...}
    // {
    // }

    template <class OElementType, class... OTags>
    inline constexpr TaggedVector(TaggedVector<OElementType, OTags...> const& other) noexcept
        : m_values {(static_cast<ElementType>(other.template get<Tags>()))...}
    {
    }

    template <
            class... Params,
            class = std::enable_if_t<((std::is_convertible_v<Params, ElementType>)&&...)>>
    inline constexpr TaggedVector(Params const&... params) noexcept
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
        using namespace detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from TaggedVector");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <class QueryTag>
    inline constexpr ElementType const& get() const noexcept
    {
        using namespace detail;
        static_assert(in_tags_v<QueryTag, tags_seq>, "requested Tag absent from TaggedVector");
        return m_values[type_seq_rank_v<QueryTag, tags_seq>];
    }

    template <std::size_t N = sizeof...(Tags)>
    constexpr inline std::enable_if_t<N == 1, ElementType const&> value() const noexcept
    {
        return m_values[0];
    }

    template <class OElementType, class... OTags>
    constexpr inline TaggedVector& operator+=(TaggedVector<OElementType, OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, detail::TypeSeq<OTags...>>);
        ((m_values[type_seq_rank_v<Tags, tags_seq>] += rhs.template get<Tags>()), ...);
        return *this;
    }

    template <class OElementType, class... OTags>
    constexpr inline TaggedVector& operator-=(TaggedVector<OElementType, OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, detail::TypeSeq<OTags...>>);
        ((m_values[type_seq_rank_v<Tags, tags_seq>] -= rhs.template get<Tags>()), ...);
        return *this;
    }

    template <class OElementType, class... OTags>
    constexpr inline TaggedVector& operator*=(TaggedVector<OElementType, OTags...> const& rhs)
    {
        static_assert(type_seq_same_v<tags_seq, detail::TypeSeq<OTags...>>);
        ((m_values[type_seq_rank_v<Tags, tags_seq>] *= rhs.template get<Tags>()), ...);
        return *this;
    }
};

template <class QueryTag, class ElementType, class... Tags>
inline constexpr ElementType const& get(TaggedVector<ElementType, Tags...> const& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class ElementType, class... Tags>
inline constexpr ElementType& get(TaggedVector<ElementType, Tags...>& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class ElementType, class... Tags, class OElementType, class... OTags>
constexpr inline auto operator+(
        TaggedVector<ElementType, Tags...> const& lhs,
        TaggedVector<OElementType, OTags...> const& rhs)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType = decltype(std::declval<ElementType>() + std::declval<OElementType>());
    return TaggedVector<RElementType, Tags...>((get<Tags>(lhs) + get<Tags>(rhs))...);
}

template <class ElementType, class... Tags, class OElementType, class... OTags>
constexpr inline auto operator-(
        TaggedVector<ElementType, Tags...> const& lhs,
        TaggedVector<OElementType, OTags...> const& rhs)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType = decltype(std::declval<ElementType>() - std::declval<OElementType>());
    return TaggedVector<RElementType, Tags...>((get<Tags>(lhs) - get<Tags>(rhs))...);
}

template <class ElementType, class... Tags, class OElementType, class... OTags>
constexpr inline auto operator*(
        TaggedVector<ElementType, Tags...> const& lhs,
        TaggedVector<OElementType, OTags...> const& rhs)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType = decltype(std::declval<ElementType>() * std::declval<OElementType>());
    return TaggedVector<RElementType, Tags...>((get<Tags>(lhs) * get<Tags>(rhs))...);
}

template <class... QueryTags, class ElementType, class... Tags>
inline constexpr TaggedVector<ElementType, QueryTags...> select(
        TaggedVector<ElementType, Tags...> const& arr) noexcept
{
    return TaggedVector<ElementType, QueryTags...>(arr);
}

template <class... QueryTags, class ElementType, class... Tags>
inline constexpr TaggedVector<ElementType, QueryTags...> select(
        TaggedVector<ElementType, Tags...>&& arr) noexcept
{
    return TaggedVector<ElementType, QueryTags...>(std::move(arr));
}

template <class ElementType>
std::ostream& operator<<(std::ostream& out, TaggedVector<ElementType> const& arr)
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
