#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

#include "taggedarray.h"

template <class, class>
class TaggedTuple;

template <class... ElementTypes, class... Tags>
class TaggedTuple<detail::TypeSeq<ElementTypes...>, detail::TypeSeq<Tags...>>
{
    std::tuple<ElementTypes...> m_values;

public:
    constexpr TaggedTuple() = default;

    constexpr TaggedTuple(TaggedTuple const&) = default;

    constexpr TaggedTuple(TaggedTuple&&) = default;

    template <class... Params>
    inline constexpr TaggedTuple(Params&&... params) noexcept
        : m_values(std::forward<Params>(params)...)
    {
    }

    template <class... OTags>
    inline constexpr TaggedTuple(TaggedTuple<OTags...> const& other) noexcept
        : m_values(other.template get<OTags>()...)
    {
    }

    template <class... OTags>
    inline constexpr TaggedTuple(TaggedTuple<OTags...>&& other) noexcept
        : m_values(std::move(other.template get<OTags>())...)
    {
    }

    constexpr inline TaggedTuple& operator=(TaggedTuple const& other) = default;

    constexpr inline TaggedTuple& operator=(TaggedTuple&& other) = default;

    template <class... OElementTypes, class... OTags>
    constexpr inline TaggedTuple& operator=(
            TaggedTuple<OElementTypes..., OTags...> const& other) noexcept
    {
        m_values = other.m_values;
        return *this;
    }

    template <class... OElementTypes, class... OTags>
    constexpr inline TaggedTuple& operator=(
            TaggedTuple<OElementTypes..., OTags...>&& other) noexcept
    {
        m_values = std::move(other.m_values);
        return *this;
    }

    constexpr inline TaggedTuple& operator=(
            std::tuple_element<0, std::tuple<ElementTypes...>> const& e) noexcept
    {
        static_assert(
                sizeof...(Tags) == 1,
                "Implicit conversion is only possible for size 1 TaggedTuples");
        std::get<0>(m_values) = e;
        return *this;
    }

    constexpr inline TaggedTuple& operator=(
            std::tuple_element<0, std::tuple<ElementTypes...>>&& e) noexcept
    {
        static_assert(
                sizeof...(Tags) == 1,
                "Implicit conversion is only possible for size 1 TaggedTuples");
        std::get<0>(m_values) = std::move(e);
        return *this;
    }

    constexpr inline operator std::tuple_element<0, std::tuple<ElementTypes...>> const &()
            const noexcept
    {
        static_assert(
                sizeof...(Tags) == 1,
                "Implicit conversion is only possible for size 1 TaggedTuples");
        return std::get<0>(m_values);
    }

    constexpr inline operator std::tuple_element<0, std::tuple<ElementTypes...>>&() noexcept
    {
        static_assert(
                sizeof...(Tags) == 1,
                "Implicit conversion is only possible for size 1 TaggedTuples");
        return std::get<0>(m_values);
    }

    template <class QueryTag>
    inline constexpr auto& get() noexcept
    {
        using namespace detail;
        return std::get<RankIn<SingleType<QueryTag>, TypeSeq<Tags...>>::val>(m_values);
    }

    template <class QueryTag>
    inline constexpr auto const& get() const noexcept
    {
        using namespace detail;
        return std::get<RankIn<SingleType<QueryTag>, TypeSeq<Tags...>>::val>(m_values);
    }
};

template <class QueryTag, class ElementType, class T>
inline constexpr auto const& get(TaggedTuple<ElementType, T> const& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}

template <class QueryTag, class ElementType, class T>
inline constexpr auto& get(TaggedTuple<ElementType, T>& tuple) noexcept
{
    return tuple.template get<QueryTag>();
}
