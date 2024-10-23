// SPDX-License-Identifier: MIT
#pragma once

namespace detail {

template <
        class TypeSeqIn,
        std::size_t Start,
        std::size_t End,
        std::size_t Idx = 0,
        class TypeSeqOut = ddc::detail::TypeSeq<>>
struct TypeSeqRange
{
    static_assert(End <= ddc::type_seq_size_v<TypeSeqIn>);
    using new_result_type = std::conditional_t<
            (Start <= Idx) && (Idx < End),
            ddc::type_seq_merge_t<
                    TypeSeqOut,
                    ddc::detail::TypeSeq<ddc::type_seq_element_t<Idx, TypeSeqIn>>>,
            TypeSeqOut>;
    using type = typename TypeSeqRange<TypeSeqIn, Start, End, Idx + 1, new_result_type>::type;
};

template <class TypeSeqIn, std::size_t Start, std::size_t End, class TypeSeqOut>
struct TypeSeqRange<TypeSeqIn, Start, End, End, TypeSeqOut>
{
    using type = TypeSeqOut;
};

} // namespace detail

/// A tool to get a subset of a TypeSeq by slicing [Start:End]
template <class TypeSeqIn, std::size_t Start, std::size_t End>
using type_seq_range_t = typename detail::TypeSeqRange<TypeSeqIn, Start, End, Start>::type;
