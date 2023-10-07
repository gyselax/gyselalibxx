// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>

namespace ddc {

namespace detail {

template <class...>
struct TypeSeq;

template <class>
struct SingleType;

template <class...>
struct TypeSeqRank;

template <class QueryTag>
struct TypeSeqRank<SingleType<QueryTag>, TypeSeq<>>
{
    static constexpr bool present = false;
};

template <class QueryTag, class... TagsTail>
struct TypeSeqRank<SingleType<QueryTag>, TypeSeq<QueryTag, TagsTail...>>
{
    static constexpr bool present = true;
    static constexpr std::size_t val = 0;
};

template <class QueryTag, class TagsHead, class... TagsTail>
struct TypeSeqRank<SingleType<QueryTag>, TypeSeq<TagsHead, TagsTail...>>
{
    static constexpr bool present
            = TypeSeqRank<SingleType<QueryTag>, TypeSeq<TagsTail...>>::present;
    static constexpr std::size_t val
            = 1 + TypeSeqRank<SingleType<QueryTag>, TypeSeq<TagsTail...>>::val;
};

template <class... QueryTags, class... Tags>
struct TypeSeqRank<TypeSeq<QueryTags...>, TypeSeq<Tags...>>
{
    using ValSeq = std::index_sequence<TypeSeqRank<QueryTags, TypeSeq<Tags...>>::val...>;
};

template <std::size_t I, class TagSeq>
struct TypeSeqElement;

template <std::size_t I, class... Tags>
struct TypeSeqElement<I, TypeSeq<Tags...>>
{
    using type = std::tuple_element_t<I, std::tuple<Tags...>>;
};

/// R contains all elements in A that are not in B.
/// Remark 1: This operation preserves the order from A.
/// Remark 2: It is similar to the set difference in the set theory (R = A\B).
/// Example: A = [a, b, c], B = [z, c, y], R = [a, b]
template <class TagSeqA, class TagSeqB, class TagSeqR>
struct TypeSeqRemove;

template <class... TagsB, class... TagsR>
struct TypeSeqRemove<TypeSeq<>, TypeSeq<TagsB...>, TypeSeq<TagsR...>>
{
    using type = TypeSeq<TagsR...>;
};

template <class HeadTagsA, class... TailTagsA, class... TagsB, class... TagsR>
struct TypeSeqRemove<TypeSeq<HeadTagsA, TailTagsA...>, TypeSeq<TagsB...>, TypeSeq<TagsR...>>
    : std::conditional_t<
              TypeSeqRank<detail::SingleType<HeadTagsA>, TypeSeq<TagsB...>>::present,
              TypeSeqRemove<TypeSeq<TailTagsA...>, TypeSeq<TagsB...>, TypeSeq<TagsR...>>,
              TypeSeqRemove<TypeSeq<TailTagsA...>, TypeSeq<TagsB...>, TypeSeq<TagsR..., HeadTagsA>>>
{
};

/// R contains all elements in A and elements in B that are not in A.
/// Remark 1: This operation preserves the order from A.
/// Remark 2: It is similar to the set union in the set theory (R = AUB).
/// Example: A = [a, b, c], B = [z, c, y], R = [a, b, c, z, y]
template <class TagSeqA, class TagSeqB, class TagSeqR>
struct TypeSeqMerge;

template <class... TagsA, class... TagsR>
struct TypeSeqMerge<TypeSeq<TagsA...>, TypeSeq<>, TypeSeq<TagsR...>>
{
    using type = TypeSeq<TagsR...>;
};

template <class... TagsA, class HeadTagsB, class... TailTagsB, class... TagsR>
struct TypeSeqMerge<TypeSeq<TagsA...>, TypeSeq<HeadTagsB, TailTagsB...>, TypeSeq<TagsR...>>
    : std::conditional_t<
              TypeSeqRank<detail::SingleType<HeadTagsB>, TypeSeq<TagsA...>>::present,
              TypeSeqMerge<TypeSeq<TagsA...>, TypeSeq<TailTagsB...>, TypeSeq<TagsR...>>,
              TypeSeqMerge<TypeSeq<TagsA...>, TypeSeq<TailTagsB...>, TypeSeq<TagsR..., HeadTagsB>>>
{
};

/// A is replaced by element of C at same position than the first element of B equal to A.
/// Remark : It may not be usefull in its own, it is an helper for TypeSeqReplace
template <class TagA, class TagSeqB, class TagSeqC>
struct TypeSeqReplaceSingle;

template <class TagA>
struct TypeSeqReplaceSingle<TagA, TypeSeq<>, TypeSeq<>>
{
    using type = TagA;
};

template <class TagA, class HeadTagsB, class... TailTagsB, class HeadTagsC, class... TailTagsC>
struct TypeSeqReplaceSingle<
        TagA,
        TypeSeq<HeadTagsB, TailTagsB...>,
        TypeSeq<HeadTagsC, TailTagsC...>>
    : std::conditional_t<
              std::is_same_v<TagA, HeadTagsB>,
              TypeSeqReplaceSingle<HeadTagsC, TypeSeq<>, TypeSeq<>>,
              TypeSeqReplaceSingle<TagA, TypeSeq<TailTagsB...>, TypeSeq<TailTagsC...>>>
{
};

/// R contains all elements of A except those of B which are replaced by those of C.
/// Remark : This operation preserves the orders.
template <class TagSeqA, class TagSeqB, class TagSeqC, class TagSeqR>
struct TypeSeqReplace;

template <class... TagsB, class... TagsC, class... TagsR>
struct TypeSeqReplace<TypeSeq<>, TypeSeq<TagsB...>, TypeSeq<TagsC...>, TypeSeq<TagsR...>>
{
    using type = TypeSeq<TagsR...>;
};

template <class HeadTagsA, class... TailTagsA, class... TagsB, class... TagsC, class... TagsR>
struct TypeSeqReplace<
        TypeSeq<HeadTagsA, TailTagsA...>,
        TypeSeq<TagsB...>,
        TypeSeq<TagsC...>,
        TypeSeq<TagsR...>>
    : TypeSeqReplace<
              TypeSeq<TailTagsA...>,
              TypeSeq<TagsB...>,
              TypeSeq<TagsC...>,
              TypeSeq<TagsR...,
                      typename TypeSeqReplaceSingle<
                              HeadTagsA,
                              TypeSeq<TagsB...>,
                              TypeSeq<TagsC...>>::type>>
{
};

} // namespace detail

template <class QueryTag, class TypeSeq>
constexpr std::size_t type_seq_rank_v = std::numeric_limits<std::size_t>::max();

template <class QueryTag, class OTypeSeq>
constexpr bool in_tags_v = false;

template <class TypeSeq, class OTypeSeq>
constexpr bool type_seq_contains_v = false;

template <class TypeSeq, class B>
constexpr bool type_seq_same_v = type_seq_contains_v<TypeSeq, B>&& type_seq_contains_v<B, TypeSeq>;

template <class QueryTag, class... Tags>
constexpr bool in_tags_v<QueryTag, detail::TypeSeq<Tags...>> = detail::
        TypeSeqRank<detail::SingleType<QueryTag>, detail::TypeSeq<Tags...>>::present;

template <class... Tags, class OTypeSeq>
constexpr bool type_seq_contains_v<
        detail::TypeSeq<Tags...>,
        OTypeSeq> = ((detail::TypeSeqRank<detail::SingleType<Tags>, OTypeSeq>::present) && ...);

template <class QueryTag, class... Tags>
constexpr std::size_t type_seq_rank_v<QueryTag, detail::TypeSeq<Tags...>> = detail::
        TypeSeqRank<detail::SingleType<QueryTag>, detail::TypeSeq<Tags...>>::val;

template <std::size_t I, class TagSeq>
using type_seq_element_t = typename detail::TypeSeqElement<I, TagSeq>::type;

template <class TagSeqA, class TagSeqB>
using type_seq_remove_t = typename detail::TypeSeqRemove<TagSeqA, TagSeqB, detail::TypeSeq<>>::type;

template <class TagSeqA, class TagSeqB>
using type_seq_merge_t = typename detail::TypeSeqMerge<TagSeqA, TagSeqB, TagSeqA>::type;

template <class TagSeqA, class TagSeqB, class TagSeqC>
using type_seq_replace_t =
        typename detail::TypeSeqReplace<TagSeqA, TagSeqB, TagSeqC, detail::TypeSeq<>>::type;
} // namespace ddc
