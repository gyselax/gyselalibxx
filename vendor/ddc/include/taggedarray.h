#pragma once

#include <array>
#include <cstddef>
#include <ostream>
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

template <class... Tags>
struct TaggedArrayPrinter;

template <class QueryTag>
struct RankIn<SingleType<QueryTag>, TypeSeq<>>
{
    static constexpr bool present = false;
};

template <class QueryTag, class... TagsTail>
struct RankIn<SingleType<QueryTag>, TypeSeq<QueryTag, TagsTail...>>
{
    static constexpr bool present = true;
    static constexpr std::size_t val = 0;
};

template <class QueryTag, class TagsHead, class... TagsTail>
struct RankIn<SingleType<QueryTag>, TypeSeq<TagsHead, TagsTail...>>
{
    static constexpr bool present = RankIn<SingleType<QueryTag>, TypeSeq<TagsTail...>>::present;
    static constexpr std::size_t val = 1 + RankIn<SingleType<QueryTag>, TypeSeq<TagsTail...>>::val;
};

template <class... QueryTags, class... Tags>
struct RankIn<TypeSeq<QueryTags...>, TypeSeq<Tags...>>
{
    using ValSeq = std::index_sequence<RankIn<QueryTags, TypeSeq<Tags...>>::val...>;
};

template <class TagsHead, class TagsNext, class... TagsTail>
struct TaggedArrayPrinter<TagsHead, TagsNext, TagsTail...>
{
    template <class ElementType, class... OTags>
    static std::ostream& print_content(
            std::ostream& out,
            TaggedArray<ElementType, OTags...> const& arr)
    {
        out << arr.template get<TagsHead>();
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
inline constexpr TaggedArray<QueryTags...> select(
        TaggedArray<ElementType, Tags...> const& arr) noexcept
{
    return TaggedArray<QueryTags...>(arr);
}

template <class... QueryTags, class ElementType, class... Tags>
inline constexpr TaggedArray<QueryTags...> select(TaggedArray<ElementType, Tags...>&& arr) noexcept
{
    return TaggedArray<QueryTags...>(std::move(arr));
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
    inline constexpr TaggedArray(Params&&... params) noexcept
        : m_values {static_cast<ElementType>(std::forward<Params>(params))...}
    {
    }

    template <class OElementType, class... OTags>
    inline constexpr TaggedArray(const TaggedArray<OElementType, OTags...>& other) noexcept
        : m_values {(::get<Tags>(other))...}
    {
    }

    template <class OElementType, class... OTags>
    inline constexpr TaggedArray(TaggedArray<OElementType, OTags...>&& other) noexcept
        : m_values {std::move(::get<Tags>(other))...}
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

    constexpr inline bool operator==(const ElementType& other) const noexcept
    {
        static_assert(
                sizeof...(Tags) == 1,
                "Implicit conversion is only possible for size 1 TaggedTuples");
        return m_values[0] == other;
    }

    constexpr inline bool operator!=(const ElementType& other) const noexcept
    {
        static_assert(
                sizeof...(Tags) == 1,
                "Implicit conversion is only possible for size 1 TaggedTuples");
        return m_values[0] != other;
    }

    constexpr inline bool operator==(const TaggedArray& other) const noexcept
    {
        return m_values == other.m_values;
    }

    constexpr inline bool operator!=(const TaggedArray& other) const noexcept
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
        static_assert(
                RankIn<SingleType<QueryTag>, TypeSeq<Tags...>>::present,
                "requested Tag absent from TaggedArray");
        return m_values[RankIn<SingleType<QueryTag>, TypeSeq<Tags...>>::val];
    }

    template <class QueryTag>
    inline constexpr ElementType const& get() const noexcept
    {
        using namespace detail;
        return m_values[RankIn<SingleType<QueryTag>, TypeSeq<Tags...>>::val];
    }
};

template <class, class>
constexpr bool has_tag_v = false;

template <class QueryTag, class ElementType, class... Tags>
constexpr bool has_tag_v<QueryTag, TaggedArray<ElementType, Tags...>> = detail::
        RankIn<detail::SingleType<QueryTag>, detail::TypeSeq<Tags...>>::present;

template <class, class>
constexpr size_t tag_rank_v = -1;

template <class QueryTag, class ElementType, class... Tags>
constexpr bool tag_rank_v<QueryTag, TaggedArray<ElementType, Tags...>> = detail::
        RankIn<detail::SingleType<QueryTag>, detail::TypeSeq<Tags...>>::val;


template <class ElementType, class... Tags>
std::ostream& operator<<(std::ostream& out, TaggedArray<ElementType, Tags...> const& arr)
{
    out << "(";
    detail::TaggedArrayPrinter<Tags...>::print_content(out, arr);
    out << ")";
    return out;
}
