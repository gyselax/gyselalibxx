#pragma once

#include <array>
#include <cstddef>
#include <utility>

template <class, class...>
class TaggedArray;

namespace detail {

template <class...>
struct TypeSeq;

template <class>
struct SingleType;

template <class...>
struct RankIn;

template <class QueryTag, class... TagsTail>
struct RankIn<SingleType<QueryTag>, TypeSeq<QueryTag, TagsTail...>>
{
    static constexpr std::size_t val = 0;
};

template <class QueryTag, class TagsHead, class... TagsTail>
struct RankIn<SingleType<QueryTag>, TypeSeq<TagsHead, TagsTail...>>
{
    static constexpr std::size_t val = 1 + RankIn<SingleType<QueryTag>, TypeSeq<TagsTail...>>::val;
};

template <class... QueryTags, class... Tags>
struct RankIn<TypeSeq<QueryTags...>, TypeSeq<Tags...>>
{
    using ValSeq = std::index_sequence<RankIn<QueryTags, TypeSeq<Tags...>>::val...>;
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


template <class ElementType, class... Tags>
class TaggedArray
{
    std::array<ElementType, sizeof...(Tags)> m_values;

public:
    constexpr TaggedArray() noexcept = default;

    constexpr TaggedArray(const TaggedArray&) noexcept = default;

    constexpr TaggedArray(TaggedArray&&) noexcept = default;

    template <class... Params>
    inline constexpr TaggedArray(Params... params) noexcept
        : m_values {std::forward<ElementType>(params)...}
    {
    }

    template <class OElementType, class... OTags>
    inline constexpr TaggedArray(const TaggedArray<OElementType, OTags...>& other) noexcept
        : m_values {get<Tags>(other)...}
    {
    }

    template <class OElementType, class... OTags>
    inline constexpr TaggedArray(TaggedArray<OElementType, OTags...>&& other) noexcept
        : m_values {std::forward<OElementType>(::get<Tags>(other))...}
    {
    }

    constexpr inline TaggedArray& operator=(const TaggedArray& other) noexcept = default;

    constexpr inline TaggedArray& operator=(TaggedArray&& other) noexcept = default;

    template <class... OTags>
    constexpr inline TaggedArray& operator=(
            const TaggedArray<ElementType, OTags...>& other) noexcept
    {
        m_values = other.m_values;
        return *this;
    }

    template <class... OTags>
    constexpr inline TaggedArray& operator=(TaggedArray<ElementType, OTags...>&& other) noexcept
    {
        m_values = std::move(other.m_values);
        return *this;
    }

    constexpr inline TaggedArray& operator=(const ElementType& e) noexcept
    {
        static_assert(
                sizeof...(Tags) == 1,
                "Implicit conversion is only possible for size 1 TaggedTuples");
        m_values = e;
        return *this;
    }

    constexpr inline TaggedArray& operator=(ElementType&& e) noexcept
    {
        static_assert(
                sizeof...(Tags) == 1,
                "Implicit conversion is only possible for size 1 TaggedTuples");
        m_values = std::move(e);
        return *this;
    }

    constexpr inline operator const ElementType&() const noexcept
    {
        static_assert(
                sizeof...(Tags) == 1,
                "Implicit conversion is only possible for size 1 TaggedArrays");
        return m_values[0];
    }

    constexpr inline operator ElementType&() noexcept
    {
        static_assert(
                sizeof...(Tags) == 1,
                "Implicit conversion is only possible for size 1 TaggedArrays");
        return m_values[0];
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
        return m_values[RankIn<SingleType<QueryTag>, TypeSeq<Tags...>>::val];
    }

    template <class QueryTag>
    inline constexpr ElementType const& get() const noexcept
    {
        using namespace detail;
        return m_values[RankIn<SingleType<QueryTag>, TypeSeq<Tags...>>::val];
    }
};
