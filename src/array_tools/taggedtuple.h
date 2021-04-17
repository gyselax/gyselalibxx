#pragma once

#include <tuple>

template <class, class...>
class TaggedTuple;

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
    static constexpr size_t val = 0;
};

template <class QueryTag, class TagsHead, class... TagsTail>
struct RankIn<SingleType<QueryTag>, TypeSeq<TagsHead, TagsTail...>>
{
    static constexpr size_t val = 1 + RankIn<SingleType<QueryTag>, TypeSeq<TagsTail...>>::val;
};

template <class... QueryTags, class... Tags>
struct RankIn<TypeSeq<QueryTags...>, TypeSeq<Tags...>>
{
    using ValSeq = std::index_sequence<RankIn<QueryTags, TypeSeq<Tags...>>::val...>;
};

} // namespace detail


template <class ElementType, class... Tags>
class TaggedTuple
{
    std::array<ElementType, sizeof...(Tags)> m_values;

public:
    constexpr TaggedTuple() noexcept = default;

    constexpr TaggedTuple(const TaggedTuple&) noexcept = default;

    constexpr TaggedTuple(TaggedTuple&&) noexcept = default;

    template <class... Params>
    inline constexpr TaggedTuple(Params... params) noexcept
        : m_values {std::forward<Params>(params)...}
    {
    }

    template <class... OTags>
    inline constexpr TaggedTuple(const TaggedTuple<OTags...>& other) noexcept
        : m_values(other.template get<OTags>()...)
    {
    }

    template <class... OTags>
    inline constexpr TaggedTuple(TaggedTuple<OTags...>&& other) noexcept
        : m_values(other.template get<OTags>()...)
    {
    }

    template <class QueryTag>
    inline constexpr auto get() noexcept
    {
        using namespace detail;
        return m_values[RankIn<SingleType<QueryTag>, TypeSeq<Tags...>>::val];
    }

    template <class QueryTag>
    inline constexpr auto get() const noexcept
    {
        using namespace detail;
        return m_values[RankIn<SingleType<QueryTag>, TypeSeq<Tags...>>::val];
    }
};
