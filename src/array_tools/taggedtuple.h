#pragma once

#include <tuple>

template <class, class...>
class TaggedTuple;

namespace detail {

/** Maps types to size_t
 */
template <class...>
struct TypeSeq;

template <class...>
struct SMapImpl;

template <class QueryTag, class... TagsTail>
struct SMapImpl<QueryTag, TypeSeq<QueryTag, TagsTail...>>
{
    static inline constexpr size_t get()
    {
        return 0;
    }
};

template <class QueryTag, class TagsHead, class... TagsTail>
struct SMapImpl<QueryTag, TypeSeq<TagsHead, TagsTail...>>
{
    static inline constexpr size_t get()
    {
        return 1 + SMapImpl<QueryTag, TypeSeq<TagsTail...>>::get();
    }
};

template <class QueryTag, class... Tags>
static inline constexpr size_t rank_of()
{
    return SMapImpl<QueryTag, Tags...>::get();
}

template <class QueryTag, class... Tags>
static inline constexpr size_t rank_of()
{
    return SMapImpl<QueryTag, Tags...>::get();
}

} // namespace detail


template <class ElementType, class... Tags>
class TaggedTuple : public std::array<ElementType, sizeof...(Tags)>
{
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

    template <class QueryTag>
    inline constexpr auto get() noexcept
    {
        return this->operator[](detail::rank_of<QueryTag, Tags...>());
    }

    template <class QueryTag>
    inline constexpr auto get() const noexcept
    {
        return this->operator[](detail::rank_of<QueryTag, Tags...>());
    }
};
