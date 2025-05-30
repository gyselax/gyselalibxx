// SPDX-License-Identifier: MIT
#pragma once

namespace detail {

template <class GridDim, typename = void>
struct GetCDim
{
    using type = GridDim;
};

template <class GridDim>
struct GetCDim<
        GridDim,
        std::enable_if_t<
                std::is_same_v<
                        typename GridDim::continuous_dimension_type,
                        typename GridDim::continuous_dimension_type>,
                void>>
{
    using type = typename GridDim::continuous_dimension_type;
};

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
            std::is_same_v<typename GetCDim<HeadGrid>::type, Dim>,
            HeadGrid,
            typename FindGrid<Dim, ddc::detail::TypeSeq<Grids...>>::type>;
};

template <class Dim>
struct FindGrid<Dim, ddc::detail::TypeSeq<>>
{
    static_assert(std::is_same_v<Dim, Dim>, "Grid not found");
    using type = void;
};

template <class CoordType, class IdxRangeType>
struct FindIdxType;

template <class... Dims, class IdxRangeType>
struct FindIdxType<Coord<Dims...>, IdxRangeType>
{
    using type = Idx<typename FindGrid<Dims, ddc::to_type_seq_t<IdxRangeType>>::type...>;
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

template <class TypeSeqType, class OrderedTypeSeqType>
struct GetPermutationParity;

template <class HeadType, class... TailTypeSeq, class... TailOrderedTypeSeq>
struct GetPermutationParity<
        ddc::detail::TypeSeq<HeadType, TailTypeSeq...>,
        ddc::detail::TypeSeq<HeadType, TailOrderedTypeSeq...>>
{
private:
    static constexpr bool contains_duplicate
            = (!std::is_same_v<
                      ddc::detail::TypeSeq<HeadType, TailOrderedTypeSeq...>,
                      typename detail::GetUnique<
                              ddc::detail::TypeSeq<HeadType, TailOrderedTypeSeq...>>::type>)
              || (!std::is_same_v<
                      ddc::detail::TypeSeq<HeadType, TailTypeSeq...>,
                      typename detail::GetUnique<
                              ddc::detail::TypeSeq<HeadType, TailTypeSeq...>>::type>);

public:
    static constexpr int value
            = contains_duplicate ? 0
                                 : GetPermutationParity<
                                         ddc::detail::TypeSeq<TailTypeSeq...>,
                                         ddc::detail::TypeSeq<TailOrderedTypeSeq...>>::value;
};

template <class HeadType, class... TailTypeSeq, class OrderedHeadType, class... TailOrderedTypeSeq>
struct GetPermutationParity<
        ddc::detail::TypeSeq<HeadType, TailTypeSeq...>,
        ddc::detail::TypeSeq<OrderedHeadType, TailOrderedTypeSeq...>>
{
private:
    static constexpr bool contains_duplicate
            = (!std::is_same_v<
                      ddc::detail::TypeSeq<OrderedHeadType, TailOrderedTypeSeq...>,
                      typename detail::GetUnique<
                              ddc::detail::TypeSeq<OrderedHeadType, TailOrderedTypeSeq...>>::type>)
              || (!std::is_same_v<
                      ddc::detail::TypeSeq<HeadType, TailTypeSeq...>,
                      typename detail::GetUnique<
                              ddc::detail::TypeSeq<HeadType, TailTypeSeq...>>::type>);

public:
    static constexpr int value
            = contains_duplicate ? 0
                                 : -GetPermutationParity<
                                         ddc::type_seq_replace_t<
                                                 ddc::detail::TypeSeq<TailTypeSeq...>,
                                                 ddc::detail::TypeSeq<OrderedHeadType>,
                                                 ddc::detail::TypeSeq<HeadType>>,
                                         ddc::detail::TypeSeq<TailOrderedTypeSeq...>>::value;
};

template <>
struct GetPermutationParity<ddc::detail::TypeSeq<>, ddc::detail::TypeSeq<>>
{
    static constexpr int value = 1;
};

template <class Element, class IndexSequence>
struct TypeSeqDuplicate;

template <class Element, std::size_t... Is>
struct TypeSeqDuplicate<Element, std::index_sequence<Is...>>
{
private:
    template <size_t>
    using get_element = Element;

public:
    using type = ddc::detail::TypeSeq<get_element<Is>...>;
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

/**
 * @brief Determine if the permutation parity of a type sequence.
 * @tparam TypeSeqType The type sequence whose permutation parity is calculated.
 * @tparam OrderedTypeSeq The final order of the indices.whose permutation parity is calculated.
 * @tparam OrderedTypeSeq The final order of the indices.
 */
template <class TypeSeqType, class OrderedTypeSeq>
constexpr int type_seq_permutation_parity_v
        = detail::GetPermutationParity<TypeSeqType, OrderedTypeSeq>::value;

/**
 * @brief Create a type sequence containing the element Element, repeated n times.
 * @tparam Element The element to be placed in the type sequence.
 * @tparam n_element The number of times the element should appear in the type sequence.
 */
template <class Element, std::size_t n_elements>
using type_seq_duplicate_t =
        typename detail::TypeSeqDuplicate<Element, std::make_index_sequence<n_elements>>::type;

/**
 * @brief Find the type of an index which allows access to a Coordinate of the specified type.
 * @tparam CoordType The type of the coordinate
 * @tparam IdxRangeType The type of the index range that the index will come from.
 */
template <class CoordType, class IdxRangeType>
using find_idx_t = typename detail::FindIdxType<CoordType, IdxRangeType>::type;
