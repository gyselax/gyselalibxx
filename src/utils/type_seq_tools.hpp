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

template <class Dim, class TypeSeqGrid>
struct FindGrid;

template <class Dim, class HeadGrid, class... Grids>
struct FindGrid<Dim, ddc::detail::TypeSeq<HeadGrid, Grids...>>
{
    using type = std::conditional_t<
            std::is_same_v<typename HeadGrid::continuous_dimension_type, Dim>,
            HeadGrid,
            typename FindGrid<Dim, ddc::detail::TypeSeq<Grids...>>::type>;
};

template <class Dim>
struct FindGrid<Dim, ddc::detail::TypeSeq<>>
{
    static_assert(std::is_same_v<Dim, Dim>, "Grid not found");
    using type = void;
};

template <class... TypeSeqs>
struct TypeSeqCat;

template <class... First, class... Second, class... TailTypeSeqs>
struct TypeSeqCat<ddc::detail::TypeSeq<First...>, ddc::detail::TypeSeq<Second...>, TailTypeSeqs...>
{
    using type =
            typename TypeSeqCat<ddc::detail::TypeSeq<First..., Second...>, TailTypeSeqs...>::type;
};

template <class... First, class... Second>
struct TypeSeqCat<ddc::detail::TypeSeq<First...>, ddc::detail::TypeSeq<Second...>>
{
    using type = ddc::detail::TypeSeq<First..., Second...>;
};

template <class StartTypeSeq, class ResultTypeSeq = ddc::detail::TypeSeq<>>
struct GetUnique;

template <class HeadElem, class... TailElems, class... ResultElems>
struct GetUnique<ddc::detail::TypeSeq<HeadElem, TailElems...>, ddc::detail::TypeSeq<ResultElems...>>
{
    using type = typename GetUnique<
            ddc::detail::TypeSeq<TailElems...>,
            std::conditional_t<
                    ddc::in_tags_v<HeadElem, ddc::detail::TypeSeq<ResultElems...>>,
                    ddc::detail::TypeSeq<ResultElems...>,
                    ddc::detail::TypeSeq<ResultElems..., HeadElem>>>::type;
};

template <class ResultTypeSeq>
struct GetUnique<ddc::detail::TypeSeq<>, ResultTypeSeq>
{
    using type = ResultTypeSeq;
};

} // namespace detail

/// A tool to get a subset of a TypeSeq by slicing [Start:End]
template <class TypeSeqIn, std::size_t Start, std::size_t End>
using type_seq_range_t = typename detail::TypeSeqRange<TypeSeqIn, Start, End, Start>::type;

/**
 * @brief Concatenate type sequences into a new type sequence.
 * This is similar to type_seq_merge_t but it does not remove duplicate elements.
 * @tparam TypeSeqs The type sequences to be concatenated.
 */
template <class... TypeSeqs>
using type_seq_cat_t = typename detail::TypeSeqCat<TypeSeqs...>::type;

/// A tool to find the grid that is defined along the specified dimension (e.g. get GridX from X)
template <class Dim, class TypeSeqGrid>
using find_grid_t = typename detail::FindGrid<Dim, TypeSeqGrid>::type;

/**
 * @brief Get a TypeSeq containing all unique types from the original TypeSeq.
 * @tparam StartTypeSeq The original TypeSeq which may contain duplicate types.
 */
template <class StartTypeSeq>
using type_seq_unique_t = typename detail::GetUnique<StartTypeSeq>::type;

/**
 * @brief Determine if a type sequence only contains unique elements.
 * @tparam TypeSeqType The type sequence being examined.
 */
template <class TypeSeqType>
constexpr bool type_seq_has_unique_elements_v
        = std::is_same_v<TypeSeqType, type_seq_unique_t<TypeSeqType>>;
