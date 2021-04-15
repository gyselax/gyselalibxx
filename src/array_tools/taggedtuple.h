#pragma once

#include <tuple>

template <class, class...>
class TaggedTuple;

namespace detail {

template <class QueryTag, class ElementType, class... TagsTail>
struct SMapImpl;

}

template <class QueryTag, class ElementType, class... Tags>
inline constexpr const ElementType& get(const TaggedTuple<ElementType, Tags...>& map) noexcept
{
    return detail::SMapImpl<QueryTag, ElementType, Tags...>::call(map);
}

template <class QueryTag, class ElementType, class... Tags>
inline constexpr ElementType& get(TaggedTuple<ElementType, Tags...>& map) noexcept
{
    return detail::SMapImpl<QueryTag, ElementType, Tags...>::call(map);
}

template <class ElementType>
class TaggedTuple<ElementType>
{
};

template <class ElementType, class TagsHead, class... TagsTail>
class TaggedTuple<ElementType, TagsHead, TagsTail...> : public TaggedTuple<ElementType, TagsTail...>
{
public:
    using Parent = TaggedTuple<ElementType, TagsTail...>;

private:
    ElementType m_value;

public:
    constexpr TaggedTuple() noexcept = default;

    template <class HeadType, class... TailTypes>
    inline constexpr TaggedTuple(HeadType head, TailTypes... tail) noexcept
        : Parent(tail...)
        , m_value(std::move(head))
    {
    }

    template <class... Params>
    inline constexpr TaggedTuple(const TaggedTuple<Params...>& other) noexcept
        : Parent(other)
        , m_value(::get<TagsHead>(other))
    {
    }

    inline constexpr ElementType& get() noexcept
    {
        return m_value;
    }

    inline constexpr const ElementType& get() const noexcept
    {
        return m_value;
    }
};

namespace detail {

template <class QueryTag, class ElementType, class... TagsTail>
struct SMapImpl<QueryTag, ElementType, QueryTag, TagsTail...>
{
    static inline constexpr const ElementType& call(
            const TaggedTuple<ElementType, QueryTag, TagsTail...>& map)
    {
        return map.get();
    }
    static inline constexpr ElementType& call(TaggedTuple<ElementType, QueryTag, TagsTail...>& map)
    {
        return map.get();
    }
};

template <class QueryTag, class ElementType, class TagsHead, class... TagsTail>
struct SMapImpl<QueryTag, ElementType, TagsHead, TagsTail...>
{
    static inline constexpr const ElementType& call(
            const TaggedTuple<ElementType, TagsHead, TagsTail...>& map)
    {
        return SMapImpl<QueryTag, ElementType, TagsTail...>::call(map);
    }
    static inline constexpr ElementType& call(TaggedTuple<ElementType, TagsHead, TagsTail...>& map)
    {
        return SMapImpl<QueryTag, ElementType, TagsTail...>::call(map);
    }
};

} // namespace detail
