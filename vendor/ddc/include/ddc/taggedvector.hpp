#pragma once

#include <array>
#include <cstddef>
#include <ostream>
#include <utility>

#include "ddc/type_seq.hpp"

template <class, class...>
class TaggedVector;

namespace detail {

template <class... Tags>
struct TaggedVectorPrinter;

template <class TagsHead, class TagsNext, class... TagsTail>
struct TaggedVectorPrinter<TagsHead, TagsNext, TagsTail...>
{
    template <class ElementType, class... OTags>
    static std::ostream& print_content(
            std::ostream& out,
            TaggedVector<ElementType, OTags...> const& arr)
    {
        out << arr.template get<TagsHead>() << ", ";
        return TaggedVectorPrinter<TagsNext, TagsTail...>::print_content(out, arr);
    }
};

template <class Tag>
struct TaggedVectorPrinter<Tag>
{
    template <class ElementType, class... OTags>
    static std::ostream& print_content(
            std::ostream& out,
            TaggedVector<ElementType, OTags...> const& arr)
    {
        out << arr.template get<Tag>();
        return out;
    }
};

template <>
struct TaggedVectorPrinter<>
{
    template <class ElementType, class... OTags>
    static std::ostream& print_content(
            std::ostream& out,
            TaggedVector<ElementType, OTags...> const& arr)
    {
        return out;
    }
};

template <class... C>
static inline constexpr void force_eval(C&&...)
{
}

template <class>
class TaggedVectorImpl;

template <class ElementType, class... Tags>
class TaggedVectorImpl<TaggedVector<ElementType, Tags...>>
{
protected:
    std::array<ElementType, sizeof...(Tags)> m_values;

public:
    static constexpr std::size_t size() noexcept
    {
        return sizeof...(Tags);
    }

public:
    inline constexpr TaggedVectorImpl() = default;

    inline constexpr TaggedVectorImpl(TaggedVectorImpl const&) = default;

    inline constexpr TaggedVectorImpl(TaggedVectorImpl&&) = default;

    template <class OElementType, class... OTags>
    inline constexpr TaggedVectorImpl(TaggedVector<OElementType, OTags...> const& other) noexcept
        : m_values {(static_cast<ElementType>(other.template get<Tags>()))...}
    {
    }

    template <class OElementType, class... OTags>
    inline constexpr TaggedVectorImpl(TaggedVector<OElementType, OTags...>&& other) noexcept
        : m_values {std::move(other.template get<Tags>())...}
    {
    }

    template <
            class... Params,
            typename ::std::enable_if<((std::is_convertible_v<Params, ElementType>)&&...), int>::
                    type
            = 0>
    inline constexpr TaggedVectorImpl(Params&&... params) noexcept
        : m_values {static_cast<ElementType>(std::forward<Params>(params))...}
    {
    }

    constexpr inline TaggedVectorImpl& operator=(TaggedVectorImpl const& other) = default;

    constexpr inline TaggedVectorImpl& operator=(TaggedVectorImpl&& other) = default;

    template <class... OTags>
    constexpr inline TaggedVectorImpl& operator=(
            TaggedVector<ElementType, OTags...> const& other) noexcept
    {
        m_values = other.m_values;
        return *this;
    }

    template <class... OTags>
    constexpr inline TaggedVectorImpl& operator=(
            TaggedVector<ElementType, OTags...>&& other) noexcept
    {
        m_values = std::move(other.m_values);
        return *this;
    }

    constexpr inline bool operator==(TaggedVectorImpl const& other) const noexcept
    {
        return m_values == other.m_values;
    }

    constexpr inline bool operator!=(TaggedVectorImpl const& other) const noexcept
    {
        return m_values != other.m_values;
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

    template <class QueryTag>
    inline constexpr ElementType& get() noexcept
    {
        using namespace detail;
        static_assert(
                in_tags_v<QueryTag, TypeSeq<Tags...>>,
                "requested Tag absent from TaggedVectorImpl");
        return m_values[type_seq_rank_v<QueryTag, TypeSeq<Tags...>>];
    }

    template <class QueryTag>
    inline constexpr ElementType const& get() const noexcept
    {
        using namespace detail;
        static_assert(
                in_tags_v<QueryTag, TypeSeq<Tags...>>,
                "requested Tag absent from TaggedVectorImpl");
        return m_values[type_seq_rank_v<QueryTag, TypeSeq<Tags...>>];
    }
};

template <class>
class SingleTagArrayImpl;

template <class ElementType, class Tag>
class SingleTagArrayImpl<TaggedVector<ElementType, Tag>>
    : public TaggedVectorImpl<TaggedVector<ElementType, Tag>>
{
public:
    inline TaggedVector<ElementType, Tag>& operator=(ElementType const& e) noexcept
    {
        this->m_values = e;
        return *this;
    }

    inline TaggedVector<ElementType, Tag>& operator=(ElementType&& e) noexcept
    {
        this->m_values = std::move(e);
        return *this;
    }

    constexpr inline bool operator==(ElementType const& other) const noexcept
    {
        return this->m_values[0] == other;
    }

    constexpr inline bool operator!=(ElementType const& other) const noexcept
    {
        return this->m_values[0] != other;
    }


    constexpr inline operator ElementType const&() const noexcept
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
inline constexpr ElementType const& get(TaggedVector<ElementType, Tags...> const& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class ElementType, class... Tags>
inline constexpr ElementType& get(TaggedVector<ElementType, Tags...>& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

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

template <class ElementType, class Tag0, class Tag1, class... Tags>
class TaggedVector<ElementType, Tag0, Tag1, Tags...>
    : public detail::TaggedVectorImpl<TaggedVector<ElementType, Tag0, Tag1, Tags...>>
{
    using Super = detail::TaggedVectorImpl<TaggedVector<ElementType, Tag0, Tag1, Tags...>>;

public:
    inline constexpr TaggedVector() = default;

    inline constexpr TaggedVector(TaggedVector const&) = default;

    inline constexpr TaggedVector(TaggedVector&&) = default;

    template <class OElementType, class... OTags>
    inline constexpr TaggedVector(TaggedVector<OElementType, OTags...> const& other) noexcept
        : Super {::get<Tag0>(other), ::get<Tag1>(other), ::get<Tags>(other)...}
    {
    }

    template <class OElementType, class... OTags>
    inline constexpr TaggedVector(TaggedVector<OElementType, OTags...>&& other) noexcept
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
    inline constexpr TaggedVector(Params&&... params) noexcept
        : Super {static_cast<ElementType>(std::forward<Params>(params))...}
    {
    }

    constexpr inline TaggedVector& operator=(TaggedVector const& other) = default;

    constexpr inline TaggedVector& operator=(TaggedVector&& other) = default;

    template <class... OTags>
    constexpr inline TaggedVector& operator=(
            TaggedVector<ElementType, OTags...> const& other) noexcept
    {
        return Super::operator=(other);
    }

    template <class... OTags>
    constexpr inline TaggedVector& operator=(TaggedVector<ElementType, OTags...>&& other) noexcept
    {
        return Super::operator=(std::move(other));
    }
};

template <class ElementType, class Tag>
class TaggedVector<ElementType, Tag>
    : public detail::SingleTagArrayImpl<TaggedVector<ElementType, Tag>>
{
    using Super = detail::SingleTagArrayImpl<TaggedVector<ElementType, Tag>>;

public:
    inline constexpr TaggedVector() = default;

    inline constexpr TaggedVector(TaggedVector const&) = default;

    inline constexpr TaggedVector(TaggedVector&&) = default;

    template <
            class OElementType,
            class... OTags,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedVector(TaggedVector<OElementType, OTags...> const& other) noexcept
        : Super {(::get<Tag>(other))}
    {
    }

    template <
            class OElementType,
            class... OTags,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedVector(TaggedVector<OElementType, OTags...>&& other) noexcept
        : Super {std::move(::get<Tag>(other))}
    {
    }

    template <
            class OElementType,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedVector(OElementType&& param) noexcept
        : Super {static_cast<ElementType>(std::forward<OElementType>(param))}
    {
    }

    constexpr inline TaggedVector& operator=(TaggedVector const& other) = default;

    constexpr inline TaggedVector& operator=(TaggedVector&& other) = default;

    template <class... OTags>
    constexpr inline TaggedVector& operator=(
            TaggedVector<ElementType, OTags...> const& other) noexcept
    {
        return Super::operator=(other);
    }

    template <class... OTags>
    constexpr inline TaggedVector& operator=(TaggedVector<ElementType, OTags...>&& other) noexcept
    {
        return Super::operator=(std::move(other));
    }
};

template <class ElementType>
class TaggedVector<ElementType> : public detail::TaggedVectorImpl<TaggedVector<ElementType>>
{
    using Super = detail::SingleTagArrayImpl<TaggedVector<ElementType>>;

public:
    inline constexpr TaggedVector() = default;

    inline constexpr TaggedVector(TaggedVector const&) = default;

    inline constexpr TaggedVector(TaggedVector&&) = default;

    template <
            class OElementType,
            class... OTags,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedVector(TaggedVector<OElementType, OTags...> const&) noexcept
    {
    }

    template <
            class OElementType,
            class... OTags,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedVector(TaggedVector<OElementType, OTags...>&&) noexcept
    {
    }

    template <
            class OElementType,
            typename ::std::enable_if_t<std::is_convertible_v<OElementType, ElementType>, int> = 0>
    inline constexpr TaggedVector(OElementType&& param) noexcept
    {
    }

    constexpr inline TaggedVector& operator=(TaggedVector const& other) = default;

    constexpr inline TaggedVector& operator=(TaggedVector&& other) = default;

    template <class... OTags>
    constexpr inline TaggedVector& operator=(
            TaggedVector<ElementType, OTags...> const& other) noexcept
    {
        return Super::operator=(other);
    }

    template <class... OTags>
    constexpr inline TaggedVector& operator=(TaggedVector<ElementType, OTags...>&& other) noexcept
    {
        return Super::operator=(std::move(other));
    }
};

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline TaggedVector<ElementType, Tags...>& operator+=(
        TaggedVector<ElementType, Tags...>& self,
        TaggedVector<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    detail::force_eval(self.template get<Tags>() += other.template get<Tags>()...);
    return self;
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline auto operator+(
        TaggedVector<ElementType, Tags...> const& one,
        TaggedVector<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType
            = decltype(std::declval<ElementType const>() + std::declval<OElementType const>());
    return TaggedVector<RElementType, Tags...>(get<Tags>(one) + get<Tags>(other)...);
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline TaggedVector<ElementType, Tags...>& operator-=(
        TaggedVector<ElementType, Tags...>& self,
        TaggedVector<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    detail::force_eval(self.template get<Tags>() -= other.template get<Tags>()...);
    return self;
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline auto operator-(
        TaggedVector<ElementType, Tags...> const& one,
        TaggedVector<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType
            = decltype(std::declval<ElementType const>() - std::declval<OElementType const>());
    return TaggedVector<RElementType, Tags...>(get<Tags>(one) - get<Tags>(other)...);
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline TaggedVector<ElementType, Tags...>& operator*=(
        TaggedVector<ElementType, Tags...>& self,
        TaggedVector<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    detail::force_eval(self.template get<Tags>() *= other.template get<Tags>()...);
    return self;
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline auto operator*(
        TaggedVector<ElementType, Tags...> const& one,
        TaggedVector<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType
            = decltype(std::declval<ElementType const>() * std::declval<OElementType const>());
    return TaggedVector<RElementType, Tags...>(get<Tags>(one) * get<Tags>(other)...);
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline TaggedVector<ElementType, Tags...>& operator/=(
        TaggedVector<ElementType, Tags...>& self,
        TaggedVector<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    detail::force_eval(self.template get<Tags>() /= other.template get<Tags>()...);
    return self;
}

template <class ElementType, class OElementType, class... Tags, class... OTags>
constexpr inline auto operator/(
        TaggedVector<ElementType, Tags...> const& one,
        TaggedVector<OElementType, OTags...> const& other)
{
    static_assert(type_seq_same_v<detail::TypeSeq<Tags...>, detail::TypeSeq<OTags...>>);
    using RElementType
            = decltype(std::declval<ElementType const>() / std::declval<OElementType const>());
    return TaggedVector<RElementType, Tags...>(get<Tags>(one) / get<Tags>(other)...);
}


template <class ElementType, class... Tags>
std::ostream& operator<<(std::ostream& out, TaggedVector<ElementType, Tags...> const& arr)
{
    out << "(";
    detail::TaggedVectorPrinter<Tags...>::print_content(out, arr);
    out << ")";
    return out;
}
