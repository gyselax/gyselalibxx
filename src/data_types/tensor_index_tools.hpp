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

template <typename Type>
inline constexpr bool is_tensor_index_element_v
        = details::enable_tensor_index_element<std::remove_const_t<std::remove_reference_t<Type>>>;

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

template <class ValidatingTensorIndexSet, class TypeSeqValidIndexSet>
struct ToTensorIndexElement;

template <class ValidatingTensorIndexSet, class... ValidIndexSet>
struct ToTensorIndexElement<ValidatingTensorIndexSet, ddc::detail::TypeSeq<ValidIndexSet...>>
{
    using type = TensorIndexElement<ValidatingTensorIndexSet, ValidIndexSet...>;
};

} // namespace details

template <class ValidatingTensorIndexSet, class TypeSeqValidIndexSet>
using to_tensor_index_element_t = typename details::
        ToTensorIndexElement<ValidatingTensorIndexSet, TypeSeqValidIndexSet>::type;

namespace details {

template <std::size_t Elem, std::size_t IdxDimHint, class TensorIndexSetType>
struct GetNthTensorIndexElement
{
    static constexpr std::size_t IdxDim = IdxDimHint - 1;
    using VectorIndexSetAlongDim = ddc::type_seq_element_t<IdxDim, TensorIndexSetType>;
    using tensor_index_type_seq = type_seq_cat_t<
            typename GetNthTensorIndexElement<
                    Elem / ddc::type_seq_size_v<VectorIndexSetAlongDim>,
                    IdxDimHint - 1,
                    TensorIndexSetType>::tensor_index_type_seq,
            ddc::detail::TypeSeq<ddc::type_seq_element_t<
                    Elem % ddc::type_seq_size_v<VectorIndexSetAlongDim>,
                    VectorIndexSetAlongDim>>>;
    using type = to_tensor_index_element_t<TensorIndexSetType, tensor_index_type_seq>;
};

template <std::size_t Elem, class TensorIndexSetType>
struct GetNthTensorIndexElement<Elem, 0, TensorIndexSetType>
{
    using tensor_index_type_seq = ddc::detail::TypeSeq<>;
};

} // namespace details

template <std::size_t Elem, class TensorIndexSetType>
using get_nth_tensor_index_element_t = typename details::GetNthTensorIndexElement<
        Elem,
        ddc::type_seq_size_v<TensorIndexSetType>,
        TensorIndexSetType>::type;

} // namespace tensor_tools
