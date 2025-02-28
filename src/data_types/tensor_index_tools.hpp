// SPDX-License-Identifier: MIT
#pragma once
#include "type_seq_tools.hpp"
#include "vector_index_tools.hpp"

namespace tensor_tools {

/**
 * @brief A class describing an index of a tensor.
 * For example for a 2x2 metric tensor on an (x,y) plane the element @f$ g_{xx} @f$
 * would have the index TensorIndexElement<TensorIndexSetXY, X, X>.
 * @tparam ValidatingTensorIndexSet The TypeSeq of VectorIndexSets in which each tag
 *          comprising this element (along each dimension) can be found.
 * @tparam Dims The dimensions at which the tensor is accessed.
 */
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
    /**
     * @brief The rank of the tensor.
     * This is equivalent to the number of indices required to
     * access an element of the tensor.
     *
     * @return The rank of the tensor.
     */
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
    /**
     * @brief Get an integer which describes the position of this element in
     * a 1D array which stores a tensor described by the TensorIndexSet.
     * @return The index of this element in a 1D array.
     */
    KOKKOS_FUNCTION static constexpr std::size_t index()
    {
        return internal_index(std::make_index_sequence<rank()>());
    }

    /**
     * @brief Get the types that describes the index in dimension IDim.
     * @tparam IDim The index of interest.
     */
    template <std::size_t IDim>
    using index_on_dim_t = ddc::type_seq_element_t<IDim, IdxTypeSeq>;

    /**
     * @brief Get the TensorIndexSet that this element indexes.
     */
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

/**
 * @brief A helper to check if a type is a TensorIndexElement.
 */
template <typename Type>
inline constexpr bool is_tensor_index_element_v
        = details::enable_tensor_index_element<std::remove_const_t<std::remove_reference_t<Type>>>;

/**
 * @brief Count the number of instances of a character in a TypeSeq of literals.
 * @tparam search_char The character that you are searching for.
 * @tparam TupleType The type of the TypeSeq of integral_constants of characters.
 */
template <char search_char, class CharTypeSeq>
constexpr std::size_t char_occurences_v = details::CountChar<search_char, CharTypeSeq>::value;

/**
 * @brief Get a TensorIndexElement from a TypeSeq of valid VectorIndexSets and a TypeSeq of indices.
 * @tparam ValidatingTensorIndexSet A TypeSeq containing the VectorIndexSets describing the tags
 *                  that can be used as indices in each dimension.
 * @tparam TypeSeqTensorIndexTag A TypeSeq containing the tags used to index the tensor.
 */
template <class ValidatingTensorIndexSet, class TypeSeqTensorIndexTag>
using to_tensor_index_element_t = typename details::
        ToTensorIndexElement<ValidatingTensorIndexSet, TypeSeqTensorIndexTag>::type;

/**
 * @brief Get a TypeSeq of valid VectorIndexSets from a TypeSeq of VectorIndexIdMaps.
 * @tparam TypeSeqVectorIndexIdMap A TypeSeq containing a VectorIndexIdMap for each dimension
 *              of the tensor.
 */
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
                    char_occurences_v<current_char, IdsFound> == 1,
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
            char_occurences_v<current_char, IdsFound> == 1,
            ddc::detail::TypeSeq<CurrentIndex, OutIndices...>,
            ddc::detail::TypeSeq<OutIndices...>>;
};

} // namespace details

//-------------------------------------------------------------------------------------------

/**
 * @brief Create a TypeSeq of VectorIndexIdMap in which each character id only appears once
 * from a TypeSeq of VectorIndexIdMaps with repeat character ids.
 * @tparam TypeSeqVectorIndexIdMap A TypeSeq containing a VectorIndexIdMap for each dimension
 *              of the tensor.
 */
template <class TypeSeqVectorIndexIdMap>
using unique_indices_t = typename details::GetUniqueIndices<TypeSeqVectorIndexIdMap>::type;

/**
 * @brief Extract the VectorIndexIdMaps whose character id only appears once in a TypeSeq of
 * VectorIndexIdMaps.
 * @tparam TypeSeqVectorIndexIdMap A TypeSeq containing a VectorIndexIdMap for each dimension
 *              of the tensor.
 */
template <class TypeSeqVectorIndexIdMap>
using non_repeated_indices_t = typename details::GetNonRepeatedIndices<
        TypeSeqVectorIndexIdMap,
        ddc::type_seq_size_v<TypeSeqVectorIndexIdMap> - 1,
        typename details::GetIndexIds<TypeSeqVectorIndexIdMap>::type>::type;

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
    using RelevantElem = typename StartElement::index_on_dim_t<RelevantElemIdx>;
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
            typename GlobalTensorIndexElement::index_on_dim_t<
                    ddc::type_seq_rank_v<ValidIndex, TypeSeqVectorIndexIdMapGlobal>>...>;
};

} // namespace details

//-------------------------------------------------------------------------------------------

/**
 * @brief Get the TensorIndexElement which indexes a Tensor at the n-th position of its internal
 * array. E.g. for a 2x2 Tensor, get_nth_tensor_index_element_t<1, TypeSeqVectorIndexSet> returns
 * the TensorIndexElement which indexes element 1 of the array, so the element {0,1} of the tensor.
 * @tparam IndexPosition The index of the underlying array for which we want to collect the
 *                      TensorIndexElement.
 * @tparam TypeSeqVectorIndexSet A TypeSeq containing the VectorIndexSets describing the valid
 *                      indices along each dimension of the tensor.
 */
template <std::size_t IndexPosition, class TypeSeqVectorIndexSet>
using get_nth_tensor_index_element_t = typename details::GetNthTensorIndexElement<
        IndexPosition,
        ddc::type_seq_size_v<TypeSeqVectorIndexSet>,
        TypeSeqVectorIndexSet>::type;

/**
 * @brief Get the n-th valid index for a tensor which is accessed according to the pattern
 * described by a TypeSeq of VectorIndexIdMaps.
 * E.g. for a 2D tensor with components A_{xx}, A_{xy}, A_{yx}, A_{yy}, indexed with
 * @code
 * TypeSeq<VectorIndexIdMap<'i', VectorIndexSet<X, Y>>, VectorIndexIdMap<'i', VectorIndexSet<X, Y>>>
 * @endcode
 *
 * - the 1st element is TensorIndexElement<TypeSeq<VectorIndexSet<X, Y>, VectorIndexSet<X, Y>>, X, X>
 * - the 2nd element is TensorIndexElement<TypeSeq<VectorIndexSet<X, Y>, VectorIndexSet<X, Y>>, Y, Y>
 *
 * A_{xy} and A_{yx} are not valid components as they do not respect the index pattern.
 *
 * @tparam Elem The element of interest.
 * @tparam TypeSeqVectorIndexIdMap A TypeSeq containing a VectorIndexIdMap for each dimension
 *              of the tensor.
 */
template <std::size_t Elem, class TypeSeqVectorIndexIdMap>
using get_nth_tensor_index_element_from_map_t = to_tensor_index_element_t<
        get_type_seq_vector_index_set_t<TypeSeqVectorIndexIdMap>,
        typename details::GetNthTensorIndexElementFromMap<
                TypeSeqVectorIndexIdMap,
                get_nth_tensor_index_element_t<
                        Elem,
                        get_type_seq_vector_index_set_t<unique_indices_t<TypeSeqVectorIndexIdMap>>>,
                ddc::type_seq_size_v<TypeSeqVectorIndexIdMap>>::type>;

/**
 * @brief Extract the relevant elements of a TensorIndexElement to create a sub-TensorIndexElement
 * using a global and a local TypeSeq of VectorIndexIdMaps to identify the relevant elements.
 * For example:
 * for GlobalTensorIndexElement = TensorIndexElement<X,Y>
 * with TypeSeqVectorIndexIdMapGlobal = TypeSeq<VectorIndexIdMap<'i', VectorIndexSet<X, Y>>, VectorIndexIdMap<'j', VectorIndexSet<X, Y>>>
 * and TypeSeqVectorIndexIdMapLocal = TypeSeq<VectorIndexIdMap<'j', VectorIndexSet<X, Y>>>
 * we obtain TensorIndexElement<Y>
 *
 * @tparam TypeSeqVectorIndexIdMapGlobal The global TypeSeq of VectorIndexIdMaps describing how the
 *              TensorIndexElement was indexed.
 * @tparam TypeSeqVectorIndexIdMapLocal The local TypeSeq of VectorIndexIdMaps describing how to
 *              identify the relevant indices for the output.
 * @tparam GlobalTensorIndexElement The starting TensorIndexElement.
 */
template <
        class TypeSeqVectorIndexIdMapGlobal,
        class TypeSeqVectorIndexIdMapLocal,
        class GlobalTensorIndexElement>
using extract_sub_tensor_element_t = typename details::ExtractSubTensorElement<
        TypeSeqVectorIndexIdMapGlobal,
        TypeSeqVectorIndexIdMapLocal,
        GlobalTensorIndexElement>::type;

} // namespace tensor_tools
