

# File tensor\_index\_tools.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**tensor\_index\_tools.hpp**](tensor__index__tools_8hpp.md)

[Go to the documentation of this file](tensor__index__tools_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include "type_seq_tools.hpp"
#include "vector_index_tools.hpp"

namespace tensor_tools {

template <class ValidatingTensorIndexSet, class... Dims>
class TensorIndexElement;

namespace details {
template <class T>
inline constexpr bool enable_tensor_index_element = false;

template <class ValidatingTensorIndexSet, class... Dims>
inline constexpr bool
        enable_tensor_index_element<TensorIndexElement<ValidatingTensorIndexSet, Dims...>> = true;
} // namespace details

template <class... ValidatingVectorIndexSets, class... Dims>
class TensorIndexElement<ddc::detail::TypeSeq<ValidatingVectorIndexSets...>, Dims...>
{
    using ValidatingTensorIndexSet = ddc::detail::TypeSeq<ValidatingVectorIndexSets...>;
    static_assert((is_vector_index_set_v<ValidatingVectorIndexSets> && ...));
    using IdxTypeSeq = ddc::detail::TypeSeq<Dims...>;

public:
    KOKKOS_FUNCTION static constexpr std::size_t rank()
    {
        return sizeof...(Dims);
    }

private:
    static_assert(
            rank() == sizeof...(ValidatingVectorIndexSets),
            "Wrong number of indices provided");

    template <std::size_t... Is>
    KOKKOS_FUNCTION static constexpr bool valid_indices(std::index_sequence<Is...>)
    {
        return ((ddc::in_tags_v<
                 ddc::type_seq_element_t<Is, IdxTypeSeq>,
                 ddc::type_seq_element_t<Is, ValidatingTensorIndexSet>>)&&...);
    }

    static_assert(
            valid_indices(std::make_index_sequence<sizeof...(Dims)>()),
            "Index is not compatible with tensor type");

    template <std::size_t... Is>
    KOKKOS_FUNCTION static constexpr std::size_t internal_index(std::index_sequence<Is...>)
    {
        Kokkos::layout_right::mapping<
                Kokkos::extents<int, ddc::type_seq_size_v<ValidatingVectorIndexSets>...>>
                mapping;
        return mapping(ddc::type_seq_rank_v<
                       ddc::type_seq_element_t<Is, IdxTypeSeq>,
                       ValidatingVectorIndexSets>...);
    }

public:
    KOKKOS_FUNCTION static constexpr std::size_t index()
    {
        return internal_index(std::make_index_sequence<rank()>());
    }

    template <std::size_t IDim>
    using index_on_dim_t = ddc::type_seq_element_t<IDim, IdxTypeSeq>;

    using tensor_index_set = ValidatingTensorIndexSet;
};

namespace details {

template <class ValidatingTensorIndexSet, class TypeSeqTensorIndexTag>
struct ToTensorIndexElement;

template <class ValidatingTensorIndexSet, class... ValidIndexSet>
struct ToTensorIndexElement<ValidatingTensorIndexSet, ddc::detail::TypeSeq<ValidIndexSet...>>
{
    using type = TensorIndexElement<ValidatingTensorIndexSet, ValidIndexSet...>;
};

template <class TypeSeqVectorIndexIdMap>
struct GetIndexIds;

template <class... ValidIndexMap>
struct GetIndexIds<ddc::detail::TypeSeq<ValidIndexMap...>>
{
    using type = ddc::detail::TypeSeq<std::integral_constant<char, ValidIndexMap::id>...>;
};

template <class TypeSeqVectorIndexIdMap>
struct ExtractTypeSeqIndexSet;

template <class... VectorIndexIdMap>
struct ExtractTypeSeqIndexSet<ddc::detail::TypeSeq<VectorIndexIdMap...>>
{
    using type = ddc::detail::TypeSeq<typename VectorIndexIdMap::possible_idx_values...>;
};

template <class TypeSeqVectorIndexMap>
struct CalculateSize;

template <class... VectorIndexIdMap>
struct CalculateSize<ddc::detail::TypeSeq<VectorIndexIdMap...>>
{
    static constexpr std::size_t value
            = (ddc::type_seq_size_v<typename VectorIndexIdMap::possible_idx_values> * ...);
};

template <char search_char, class CharTypeSeq>
struct CountChar;

template <char search_char>
struct CountChar<search_char, ddc::detail::TypeSeq<>>
{
    static constexpr std::size_t value = 0;
};

template <char search_char, class head_char_lit, class... char_lit>
struct CountChar<search_char, ddc::detail::TypeSeq<head_char_lit, char_lit...>>
{
    static constexpr std::size_t value
            = (head_char_lit::value == search_char)
              + CountChar<search_char, ddc::detail::TypeSeq<char_lit...>>::value;
};

} // namespace details

//-------------------------------------------------------------------------------------------

template <typename Type>
inline constexpr bool is_tensor_index_element_v
        = details::enable_tensor_index_element<std::remove_const_t<std::remove_reference_t<Type>>>;

template <char search_char, class CharTypeSeq>
constexpr std::size_t char_occurrences_v = details::CountChar<search_char, CharTypeSeq>::value;

template <class ValidatingTensorIndexSet, class TypeSeqTensorIndexTag>
using to_tensor_index_element_t = typename details::
        ToTensorIndexElement<ValidatingTensorIndexSet, TypeSeqTensorIndexTag>::type;

template <class TypeSeqVectorIndexIdMap>
using get_type_seq_vector_index_set_t =
        typename details::ExtractTypeSeqIndexSet<TypeSeqVectorIndexIdMap>::type;

//-------------------------------------------------------------------------------------------

namespace details {

template <class StartTypeSeq, class ResultTypeSeq = ddc::detail::TypeSeq<>>
struct GetUniqueIndices;

template <class HeadElem, class... TailElems, class... ResultElems>
struct GetUniqueIndices<
        ddc::detail::TypeSeq<HeadElem, TailElems...>,
        ddc::detail::TypeSeq<ResultElems...>>
{
    using InsertionType = VectorIndexIdMap<
            HeadElem::id,
            get_contravariant_dims_t<typename HeadElem::possible_idx_values>>;
    using type = typename GetUniqueIndices<
            ddc::detail::TypeSeq<TailElems...>,
            std::conditional_t<
                    ddc::in_tags_v<InsertionType, ddc::detail::TypeSeq<ResultElems...>>,
                    ddc::detail::TypeSeq<ResultElems...>,
                    ddc::detail::TypeSeq<ResultElems..., InsertionType>>>::type;
};

template <class ResultTypeSeq>
struct GetUniqueIndices<ddc::detail::TypeSeq<>, ResultTypeSeq>
{
    using type = ResultTypeSeq;
};

template <
        class TypeSeqVectorIndexIdMap,
        std::size_t Elem,
        class IdsFound,
        class ResultTuple = ddc::detail::TypeSeq<>>
struct GetNonRepeatedIndices;

template <class TypeSeqVectorIndexIdMap, std::size_t Elem, class IdsFound, class... OutIndices>
struct GetNonRepeatedIndices<
        TypeSeqVectorIndexIdMap,
        Elem,
        IdsFound,
        ddc::detail::TypeSeq<OutIndices...>>
{
    using CurrentIndex = ddc::type_seq_element_t<Elem, TypeSeqVectorIndexIdMap>;
    static constexpr char current_char = CurrentIndex::id;
    static constexpr std::size_t n_count = CountChar<current_char, IdsFound>::value == 1;

    using type = typename GetNonRepeatedIndices<
            TypeSeqVectorIndexIdMap,
            Elem - 1,
            IdsFound,
            std::conditional_t<
                    char_occurrences_v<current_char, IdsFound> == 1,
                    ddc::detail::TypeSeq<CurrentIndex, OutIndices...>,
                    ddc::detail::TypeSeq<OutIndices...>>>::type;
};

template <class TypeSeqVectorIndexIdMap, class IdsFound, class... OutIndices>
struct GetNonRepeatedIndices<
        TypeSeqVectorIndexIdMap,
        0,
        IdsFound,
        ddc::detail::TypeSeq<OutIndices...>>
{
    using CurrentIndex = ddc::type_seq_element_t<0, TypeSeqVectorIndexIdMap>;
    static constexpr char current_char = CurrentIndex::id;
    using type = std::conditional_t<
            char_occurrences_v<current_char, IdsFound> == 1,
            ddc::detail::TypeSeq<CurrentIndex, OutIndices...>,
            ddc::detail::TypeSeq<OutIndices...>>;
};

template <char ID, class TypeSeqVectorIndexIdMap, class ResultTypeSeq = ddc::detail::TypeSeq<>>
struct GetRelevantVectorIndexSets;

template <
        char ID,
        class HeadVectorIndexIdMap,
        class... TailVectorIndexIdMap,
        class... OutVectorIndexSets>
struct GetRelevantVectorIndexSets<
        ID,
        ddc::detail::TypeSeq<HeadVectorIndexIdMap, TailVectorIndexIdMap...>,
        ddc::detail::TypeSeq<OutVectorIndexSets...>>
{
    using type = typename GetRelevantVectorIndexSets<
            ID,
            ddc::detail::TypeSeq<TailVectorIndexIdMap...>,
            std::conditional_t<
                    HeadVectorIndexIdMap::id == ID,
                    ddc::detail::TypeSeq<
                            typename HeadVectorIndexIdMap::possible_idx_values,
                            OutVectorIndexSets...>,
                    ddc::detail::TypeSeq<OutVectorIndexSets...>>>::type;
};

template <char ID, class ResultTypeSeq>
struct GetRelevantVectorIndexSets<ID, ddc::detail::TypeSeq<>, ResultTypeSeq>
{
    using type = ResultTypeSeq;
};

} // namespace details

//-------------------------------------------------------------------------------------------
template <class TypeSeqVectorIndexIdMap>
using index_identifiers_t = typename details::GetIndexIds<TypeSeqVectorIndexIdMap>::type;

template <class TypeSeqVectorIndexIdMap>
using unique_indices_t = typename details::GetUniqueIndices<TypeSeqVectorIndexIdMap>::type;

template <class TypeSeqVectorIndexIdMap>
using non_repeated_indices_t = typename details::GetNonRepeatedIndices<
        TypeSeqVectorIndexIdMap,
        ddc::type_seq_size_v<TypeSeqVectorIndexIdMap> - 1,
        index_identifiers_t<TypeSeqVectorIndexIdMap>>::type;

template <char ID, class TypeSeqVectorIndexIdMap>
using relevant_vector_index_sets_t =
        typename details::GetRelevantVectorIndexSets<ID, TypeSeqVectorIndexIdMap>::type;

//-------------------------------------------------------------------------------------------

namespace details {

template <std::size_t Elem, std::size_t IdxDimHint, class TypeSeqVectorIndexSet>
struct GetNthTensorIndexElement
{
    static constexpr std::size_t IdxDim = IdxDimHint - 1;
    using VectorIndexSetAlongDim = ddc::type_seq_element_t<IdxDim, TypeSeqVectorIndexSet>;
    using tensor_index_type_seq = type_seq_cat_t<
            typename GetNthTensorIndexElement<
                    Elem / ddc::type_seq_size_v<VectorIndexSetAlongDim>,
                    IdxDimHint - 1,
                    TypeSeqVectorIndexSet>::tensor_index_type_seq,
            ddc::detail::TypeSeq<ddc::type_seq_element_t<
                    Elem % ddc::type_seq_size_v<VectorIndexSetAlongDim>,
                    VectorIndexSetAlongDim>>>;
    using type = to_tensor_index_element_t<TypeSeqVectorIndexSet, tensor_index_type_seq>;
};

template <std::size_t Elem, class TypeSeqVectorIndexSet>
struct GetNthTensorIndexElement<Elem, 0, TypeSeqVectorIndexSet>
{
    using tensor_index_type_seq = ddc::detail::TypeSeq<>;
};

template <
        class TypeSeqVectorIndexIdMap,
        class StartElement,
        std::size_t Dim,
        class ResultTypeSeq = ddc::detail::TypeSeq<>>
struct GetNthTensorIndexElementFromMap;

template <class TypeSeqVectorIndexIdMap, class StartElement, std::size_t Dim, class... ResultElems>
struct GetNthTensorIndexElementFromMap<
        TypeSeqVectorIndexIdMap,
        StartElement,
        Dim,
        ddc::detail::TypeSeq<ResultElems...>>
{
    using ValidIndex = ddc::type_seq_element_t<Dim - 1, TypeSeqVectorIndexIdMap>;
    using ContraValidIndex = VectorIndexIdMap<
            ValidIndex::id,
            get_contravariant_dims_t<typename ValidIndex::possible_idx_values>>;
    static constexpr std::size_t RelevantElemIdx
            = ddc::type_seq_rank_v<ContraValidIndex, unique_indices_t<TypeSeqVectorIndexIdMap>>;
    using RelevantElem = typename StartElement::template index_on_dim_t<RelevantElemIdx>;
    using Elem = std::conditional_t<
            is_contravariant_vector_index_set_v<typename ValidIndex::possible_idx_values>,
            RelevantElem,
            typename RelevantElem::Dual>;
    using type = typename GetNthTensorIndexElementFromMap<
            TypeSeqVectorIndexIdMap,
            StartElement,
            Dim - 1,
            ddc::detail::TypeSeq<Elem, ResultElems...>>::type;
};

template <class TypeSeqVectorIndexIdMap, class StartElement, class... ResultElems>
struct GetNthTensorIndexElementFromMap<
        TypeSeqVectorIndexIdMap,
        StartElement,
        0,
        ddc::detail::TypeSeq<ResultElems...>>
{
    using type = ddc::detail::TypeSeq<ResultElems...>;
};

template <
        class TypeSeqVectorIndexIdMapGlobal,
        class TypeSeqVectorIndexIdMapLocal,
        class GlobalTensorIndexElement>
struct ExtractSubTensorElement;

template <class TypeSeqVectorIndexIdMapGlobal, class... ValidIndex, class GlobalTensorIndexElement>
struct ExtractSubTensorElement<
        TypeSeqVectorIndexIdMapGlobal,
        ddc::detail::TypeSeq<ValidIndex...>,
        GlobalTensorIndexElement>
{
    using type = TensorIndexElement<
            ddc::detail::TypeSeq<typename ValidIndex::possible_idx_values...>,
            typename GlobalTensorIndexElement::template index_on_dim_t<
                    ddc::type_seq_rank_v<ValidIndex, TypeSeqVectorIndexIdMapGlobal>>...>;
};

} // namespace details

//-------------------------------------------------------------------------------------------

template <std::size_t IndexPosition, class TypeSeqVectorIndexSet>
using get_nth_tensor_index_element_t = typename details::GetNthTensorIndexElement<
        IndexPosition,
        ddc::type_seq_size_v<TypeSeqVectorIndexSet>,
        TypeSeqVectorIndexSet>::type;

template <std::size_t Elem, class TypeSeqVectorIndexIdMap>
using get_nth_tensor_index_element_from_map_t = to_tensor_index_element_t<
        get_type_seq_vector_index_set_t<TypeSeqVectorIndexIdMap>,
        typename details::GetNthTensorIndexElementFromMap<
                TypeSeqVectorIndexIdMap,
                get_nth_tensor_index_element_t<
                        Elem,
                        get_type_seq_vector_index_set_t<unique_indices_t<TypeSeqVectorIndexIdMap>>>,
                ddc::type_seq_size_v<TypeSeqVectorIndexIdMap>>::type>;

template <
        class TypeSeqVectorIndexIdMapGlobal,
        class TypeSeqVectorIndexIdMapLocal,
        class GlobalTensorIndexElement>
using extract_sub_tensor_element_t = typename details::ExtractSubTensorElement<
        TypeSeqVectorIndexIdMapGlobal,
        TypeSeqVectorIndexIdMapLocal,
        GlobalTensorIndexElement>::type;

} // namespace tensor_tools
```


