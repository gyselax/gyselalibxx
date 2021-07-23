#pragma once

#include <array>
#include <cstddef>
#include <ostream>
#include <utility>

#include "type_seq.h"

template <class, class...>
class TaggedArray;

namespace detail {

template <class... Tags>
struct TaggedArrayPrinter;

template <class TagsHead, class TagsNext, class... TagsTail>
struct TaggedArrayPrinter<TagsHead, TagsNext, TagsTail...>
{
    template <class ElementType, class... OTags>
    static std::ostream& print_content(
            std::ostream& out,
            TaggedArray<ElementType, OTags...> const& arr)
    {
        out << arr.template get<TagsHead>() << ", ";
        return TaggedArrayPrinter<TagsNext, TagsTail...>::print_content(out, arr);
    }
};

template <class Tag>
struct TaggedArrayPrinter<Tag>
{
    template <class ElementType, class... OTags>
    static std::ostream& print_content(
            std::ostream& out,
            TaggedArray<ElementType, OTags...> const& arr)
    {
        out << arr.template get<Tag>();
        return out;
    }
};

template <>
struct TaggedArrayPrinter<>
{
    template <class ElementType, class... OTags>
    static std::ostream& print_content(
            std::ostream& out,
            TaggedArray<ElementType, OTags...> const& arr)
    {
        return out;
    }
};

template <class... C>
static inline constexpr void force_eval(C&&...)
{
}

template <class>
class TaggedArrayImpl;

template <class ElementType, class... Tags>
class TaggedArrayImpl<TaggedArray<ElementType, Tags...>>
{
protected:
    std::array<ElementType, sizeof...(Tags)> m_values;

public:
    static constexpr std::size_t size() noexcept
    {
        return sizeof...(Tags);
    }

public:
    inline constexpr TaggedArrayImpl() = default;

    inline constexpr TaggedArrayImpl(TaggedArrayImpl const&) = default;

    inline constexpr TaggedArrayImpl(TaggedArrayImpl&&) = default;

    template <class OElementType, class... OTags>
    inline constexpr TaggedArrayImpl(TaggedArray<OElementType, OTags...> const& other) noexcept
        : m_values {(other.template get<Tags>())...}
    {
    }

    template <class OElementType, class... OTags>
    inline constexpr TaggedArrayImpl(TaggedArray<OElementType, OTags...>&& other) noexcept
        : m_values {std::move(other.template get<Tags>())...}
    {
    }

    template <
            class... Params,
            typename ::std::enable_if<((std::is_convertible_v<Params, ElementType>)&&...), int>::
                    type
            = 0>
    inline constexpr TaggedArrayImpl(Params&&... params) noexcept
        : m_values {static_cast<ElementType>(std::forward<Params>(params))...}
    {
    }

    constexpr inline TaggedArrayImpl& operator=(const TaggedArrayImpl& other) = default;

    constexpr inline TaggedArrayImpl& operator=(TaggedArrayImpl&& other) = default;

    template <class... OTags>
    constexpr inline TaggedArrayImpl& operator=(
            const TaggedArray<ElementType, OTags...>& other) noexcept
    {
        m_values = other.m_values;
        return *this;
    }

    template <class... OTags>
    constexpr inline TaggedArrayImpl& operator=(TaggedArray<ElementType, OTags...>&& other) noexcept
    {
        m_values = std::move(other.m_values);
        return *this;
    }

    constexpr inline bool operator==(const TaggedArrayImpl& other) const noexcept
    {
        return m_values == other.m_values;
    }

    constexpr inline bool operator!=(const TaggedArrayImpl& other) const noexcept
    {
        return m_values != other.m_values;
    }

    /// Returns a reference to the underlying `std::array`
    constexpr inline std::array<ElementType, sizeof...(Tags)>& array() noexcept
    {
        return m_values;
    }

    /// Returns a const reference to the underlying `std::array`
    constexpr inline const std::array<ElementType, sizeof...(Tags)>& array() const noexcept
    {
        return m_values;
    }

    constexpr inline ElementType& operator[](size_t pos)
    {
        return m_values[pos];
    }

    constexpr inline const ElementType& operator[](size_t pos) const
    {
        return m_values[pos];
    }

    template <class QueryTag>
    inline constexpr ElementType& get() noexcept
    {
        using namespace detail;
        static_assert(
                in_tags_v<QueryTag, TypeSeq<Tags...>>,
                "requested Tag absent from TaggedArrayImpl");
        return m_values[type_seq_rank_v<QueryTag, TypeSeq<Tags...>>];
    }

    template <class QueryTag>
    inline constexpr ElementType const& get() const noexcept
    {
        using namespace detail;
        static_assert(
                in_tags_v<QueryTag, TypeSeq<Tags...>>,
                "requested Tag absent from TaggedArrayImpl");
        return m_values[type_seq_rank_v<QueryTag, TypeSeq<Tags...>>];
    }
};

template <class>
class SingleTagArrayImpl;

template <class ElementType, class Tag>
class SingleTagArrayImpl<TaggedArray<ElementType, Tag>>
    : public TaggedArrayImpl<TaggedArray<ElementType, Tag>>
{
public:
    inline TaggedArray<ElementType, Tag>& operator=(const ElementType& e) noexcept
    {
        this->m_values = e;
        return *this;
    }

    inline TaggedArray<ElementType, Tag>& operator=(ElementType&& e) noexcept
    {
        this->m_values = std::move(e);
        return *this;
    }

    constexpr inline bool operator==(const ElementType& other) const noexcept
    {
        return this->m_values[0] == other;
    }

    constexpr inline bool operator!=(const ElementType& other) const noexcept
    {
        return this->m_values[0] != other;
    }


    constexpr inline operator const ElementType&() const noexcept
    {
        return this->m_values[0];
    }

    constexpr inline operator ElementType&() noexcept
    {
        return this->m_values[0];
    }
};


} // namespace detail


template <class QueryTag, class ElementType, class... Tags>
inline constexpr ElementType const& get(TaggedArray<ElementType, Tags...> const& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class ElementType, class... Tags>
inline constexpr ElementType& get(TaggedArray<ElementType, Tags...>& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class... QueryTags, class ElementType, class... Tags>
inline constexpr TaggedArray<ElementType, QueryTags...> select(
        TaggedArray<ElementType, Tags...> const& arr) noexcept
{
    return TaggedArray<ElementType, QueryTags...>(arr);
}

template <class... QueryTags, class ElementType, class... Tags>
inline constexpr TaggedArray<ElementType, QueryTags...> select(
        TaggedArray<ElementType, Tags...>&& arr) noexcept
{
    return TaggedArray<ElementType, QueryTags...>(std::move(arr));
}

template <class ElementType, class Tag0, class Tag1, class... Tags>
class TaggedArray<ElementType, Tag0, Tag1, Tags...>
    : public detail::TaggedArrayImpl<TaggedArray<ElementType, Tag0, Tag1, Tags...>>
{
    using Super = detail::TaggedArrayImpl<TaggedArray<ElementType, Tag0, Tag1, Tags...>>;

public:
    inline constexpr TaggedArray() = default;

    inline constexpr TaggedArray(TaggedArray const&) = default;

    inline constexpr TaggedArray(TaggedArray&&) = default;

    template <class OElementType, class... OTags>
    inline constexpr TaggedArray(TaggedArray<OElementType, OTags...> const& other) noexcept
        : Super {::get<Tag0>(other), ::get<Tag1>(other), ::get<Tags>(other)...}
    {
    }

    template <class OElementType, class... OTags>
    inline constexpr TaggedArray(TaggedArray<OElementType, OTags...>&& other) noexcept
        : Super {
                std::move(::get<Tag0>(other)),
                std::move(::get<Tag1>(other)),
                std::move(::get<Tags>(other))...}
    {
    }

    template <
            class... Params,
            typename ::std::enable_if<((std::is_convertible_v<Params, ElementType>)&&...), int>::
                    type
            = 0>
    inline constexpr TaggedArray(Params&&... params) noexcept
        : Super {static_cast<ElementType>(std::forward<Params>(params))...}
    {
    }

    inline constexpr TaggedArray& operator=(TaggedArray const&) = default;

    inline constexpr TaggedArray& operator=(TaggedArray&&) = default;
};

template <class ElementType, class Tag>
class TaggedArray<ElementType, Tag>
    : public detail::SingleTagArrayImpl<TaggedArray<ElementType, Tag>>
{
    using Super = detail::SingleTagArrayImpl<TaggedArray<ElementType, Tag>>;

public:
    inline constexpr TaggedArray() = default;

    inline constexpr TaggedArray(TaggedArray const&) = default;

    inline constexpr TaggedArray(TaggedArray&&) = default;

    template <
            class OElementType,
            class... OTags,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedArray(TaggedArray<OElementType, OTags...> const& other) noexcept
        : Super {(::get<Tag>(other))}
    {
    }

    template <
            class OElementType,
            class... OTags,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedArray(TaggedArray<OElementType, OTags...>&& other) noexcept
        : Super {std::move(::get<Tag>(other))}
    {
    }

    template <
            class OElementType,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedArray(OElementType&& param) noexcept
        : Super {static_cast<ElementType>(std::forward<OElementType>(param))}
    {
    }
};

template <class ElementType>
class TaggedArray<ElementType> : public detail::TaggedArrayImpl<TaggedArray<ElementType>>
{
    using Super = detail::SingleTagArrayImpl<TaggedArray<ElementType>>;

public:
    inline constexpr TaggedArray() = default;

    inline constexpr TaggedArray(TaggedArray const&) = default;

    inline constexpr TaggedArray(TaggedArray&&) = default;

    template <
            class OElementType,
            class... OTags,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedArray(TaggedArray<OElementType, OTags...> const&) noexcept
    {
    }

    template <
            class OElementType,
            class... OTags,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedArray(TaggedArray<OElementType, OTags...>&&) noexcept
    {
    }

    template <
            class OElementType,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedArray(OElementType&& param) noexcept
    {
    }
};

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline TaggedArray<ElementType, Tags...>& operator+=(
        TaggedArray<ElementType, Tags...>& self,
        TaggedArray<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    detail::force_eval(self.template get<Tags>() += other.template get<Tags>()...);
    return self;
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline auto operator+(
        TaggedArray<ElementType, Tags...> const& one,
        TaggedArray<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType
            = decltype(std::declval<ElementType const>() + std::declval<OElementType const>());
    return TaggedArray<RElementType, Tags...>(get<Tags>(one) + get<Tags>(other)...);
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline TaggedArray<ElementType, Tags...>& operator-=(
        TaggedArray<ElementType, Tags...>& self,
        TaggedArray<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    detail::force_eval(self.template get<Tags>() -= other.template get<Tags>()...);
    return self;
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline auto operator-(
        TaggedArray<ElementType, Tags...> const& one,
        TaggedArray<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType
            = decltype(std::declval<ElementType const>() - std::declval<OElementType const>());
    return TaggedArray<RElementType, Tags...>(get<Tags>(one) - get<Tags>(other)...);
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline TaggedArray<ElementType, Tags...>& operator*=(
        TaggedArray<ElementType, Tags...>& self,
        TaggedArray<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    detail::force_eval(self.template get<Tags>() *= other.template get<Tags>()...);
    return self;
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline auto operator*(
        TaggedArray<ElementType, Tags...> const& one,
        TaggedArray<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType
            = decltype(std::declval<ElementType const>() * std::declval<OElementType const>());
    return TaggedArray<RElementType, Tags...>(get<Tags>(one) * get<Tags>(other)...);
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline TaggedArray<ElementType, Tags...>& operator/=(
        TaggedArray<ElementType, Tags...>& self,
        TaggedArray<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    detail::force_eval(self.template get<Tags>() /= other.template get<Tags>()...);
    return self;
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline auto operator/(
        TaggedArray<ElementType, Tags...> const& one,
        TaggedArray<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType
            = decltype(std::declval<ElementType const>() / std::declval<OElementType const>());
    return TaggedArray<RElementType, Tags...>(get<Tags>(one) / get<Tags>(other)...);
}


template <class ElementType, class... Tags>
std::ostream& operator<<(std::ostream& out, TaggedArray<ElementType, Tags...> const& arr)
{
    out << "(";
    detail::TaggedArrayPrinter<Tags...>::print_content(out, arr);
    out << ")";
    return out;
}
