// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>

namespace ddc {

namespace ddc_detail {

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
              TypeSeqRank<ddc_detail::SingleType<HeadTagsA>, TypeSeq<TagsB...>>::present,
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
              TypeSeqRank<ddc_detail::SingleType<HeadTagsB>, TypeSeq<TagsA...>>::present,
              TypeSeqMerge<TypeSeq<TagsA...>, TypeSeq<TailTagsB...>, TypeSeq<TagsR...>>,
              TypeSeqMerge<TypeSeq<TagsA...>, TypeSeq<TailTagsB...>, TypeSeq<TagsR..., HeadTagsB>>>
{
};

} // namespace ddc_detail

template <class QueryTag, class TypeSeq>
constexpr std::size_t type_seq_rank_v = std::numeric_limits<std::size_t>::max();

template <class QueryTag, class OTypeSeq>
constexpr bool in_tags_v = false;

template <class TypeSeq, class OTypeSeq>
constexpr bool type_seq_contains_v = false;

template <class TypeSeq, class B>
constexpr bool type_seq_same_v = type_seq_contains_v<TypeSeq, B>&& type_seq_contains_v<B, TypeSeq>;

template <class QueryTag, class... Tags>
constexpr bool in_tags_v<QueryTag, ddc_detail::TypeSeq<Tags...>> = ddc_detail::
        TypeSeqRank<ddc_detail::SingleType<QueryTag>, ddc_detail::TypeSeq<Tags...>>::present;

template <class... Tags, class OTypeSeq>
constexpr bool type_seq_contains_v<
        ddc_detail::TypeSeq<Tags...>,
        OTypeSeq> = ((ddc_detail::TypeSeqRank<ddc_detail::SingleType<Tags>, OTypeSeq>::present) && ...);

template <class QueryTag, class... Tags>
constexpr std::size_t type_seq_rank_v<QueryTag, ddc_detail::TypeSeq<Tags...>> = ddc_detail::
        TypeSeqRank<ddc_detail::SingleType<QueryTag>, ddc_detail::TypeSeq<Tags...>>::val;

template <std::size_t I, class TagSeq>
using type_seq_element_t = typename ddc_detail::TypeSeqElement<I, TagSeq>::type;

template <class TagSeqA, class TagSeqB>
using type_seq_remove_t =
        typename ddc_detail::TypeSeqRemove<TagSeqA, TagSeqB, ddc_detail::TypeSeq<>>::type;

template <class TagSeqA, class TagSeqB>
using type_seq_merge_t = typename ddc_detail::TypeSeqMerge<TagSeqA, TagSeqB, TagSeqA>::type;

} // namespace ddc
