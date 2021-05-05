#pragma once

#include <cstddef>
#include <tuple>
#include <utility>


template <class, class>
class TaggedTuple;

template <class...>
struct TypeSeq;

namespace detail
{

template <class>
struct SingleType;

template <class...>
struct RankIn;

template <class QueryTag, class... TagsTail>
struct RankIn<SingleType<QueryTag>, TypeSeq<QueryTag, TagsTail...>> {
    static constexpr std::size_t val = 0;
};

template <class QueryTag, class TagsHead, class... TagsTail>
struct RankIn<SingleType<QueryTag>, TypeSeq<TagsHead, TagsTail...>> {
    static constexpr std::size_t val = 1 + RankIn<SingleType<QueryTag>, TypeSeq<TagsTail...>>::val;
};

template <class... QueryTags, class... Tags>
struct RankIn<TypeSeq<QueryTags...>, TypeSeq<Tags...>> {
    using ValSeq = std::index_sequence<RankIn<QueryTags, TypeSeq<Tags...>>::val...>;
};

} // namespace detail


template <class... ElementTypes, class... Tags>
class TaggedTuple<TypeSeq<ElementTypes...>, TypeSeq<Tags...>>
{
    std::tuple<ElementTypes...> m_values;

public:
    constexpr TaggedTuple() noexcept = default;

    constexpr TaggedTuple ( const TaggedTuple& ) noexcept = default;

    constexpr TaggedTuple ( TaggedTuple&& ) noexcept = default;

    template <class... Params>
    inline constexpr TaggedTuple ( Params... params ) noexcept
        : m_values
    {
        std::forward<Params> ( params )...
    }
    {
    }

    template <class... OTags>
    inline constexpr TaggedTuple ( const TaggedTuple<OTags...>& other ) noexcept
        : m_values ( other.template get<OTags>()... )
    {
    }

    template <class... OTags>
    inline constexpr TaggedTuple ( TaggedTuple<OTags...>&& other ) noexcept
        : m_values ( other.template get<OTags>()... )
    {
    }

    constexpr inline TaggedTuple& operator= ( const TaggedTuple& other ) noexcept = default;

    constexpr inline TaggedTuple& operator= ( TaggedTuple&& other ) noexcept = default;

    template <class... OElementTypes, class... OTags>
    constexpr inline TaggedTuple& operator= (
        const TaggedTuple<OElementTypes..., OTags...>& other ) noexcept
    {
        m_values = other.m_values;
        return *this;
    }

    template <class... OElementTypes, class... OTags>
    constexpr inline TaggedTuple& operator= ( TaggedTuple<OElementTypes..., OTags...>&& other ) noexcept
    {
        m_values = std::move ( other.m_values );
        return *this;
    }

    constexpr inline TaggedTuple& operator= ( const std::tuple_element<0, std::tuple<ElementTypes...>>& e ) noexcept
    {
        static_assert (
            sizeof... ( Tags ) == 1,
            "Implicit conversion is only possible for size 1 TaggedTuples" );
        std::get<0> ( m_values ) = e;
        return *this;
    }

    constexpr inline TaggedTuple& operator= ( std::tuple_element<0, std::tuple<ElementTypes...>>&& e ) noexcept
    {
        static_assert (
            sizeof... ( Tags ) == 1,
            "Implicit conversion is only possible for size 1 TaggedTuples" );
        std::get<0> ( m_values ) = std::move ( e );
        return *this;
    }

    constexpr inline operator const std::tuple_element<0, std::tuple<ElementTypes...>>&() const noexcept
    {
        static_assert (
            sizeof... ( Tags ) == 1,
            "Implicit conversion is only possible for size 1 TaggedTuples" );
        return std::get<0> ( m_values );
    }

    constexpr inline operator std::tuple_element<0, std::tuple<ElementTypes...>>&() noexcept
    {
        static_assert (
            sizeof... ( Tags ) == 1,
            "Implicit conversion is only possible for size 1 TaggedTuples" );
        return std::get<0> ( m_values );
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

template <class QueryTag, class ElementType, class... Tags>
inline constexpr auto get ( const TaggedTuple<ElementType, Tags...>& tuple ) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class ElementType, class... Tags>
inline constexpr auto get ( TaggedTuple<ElementType, Tags...>& tuple ) noexcept
{
    return tuple.template get<QueryTag>();
}

