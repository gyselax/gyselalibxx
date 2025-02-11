// SPDX-License-Identifier: MIT
#pragma once

#include "tensor.hpp"
#include "tensor_index_set.hpp"
#include "vector_index_set.hpp"

namespace tensor_tools {

template <class T>
inline constexpr bool enable_indexed_tensor = false;

template <class TensorType, class TensorIndexIdMapType>
class IndexedTensor
{
public:
    using index_pattern = TensorIndexIdMapType;
    using tensor_type = TensorType;

private:
    TensorType& m_tensor;

public:
    IndexedTensor(TensorType& tensor) : m_tensor(tensor) {}

    TensorType& operator()()
    {
        return m_tensor;
    }
};

template <class TensorType, class TensorIdxIdMapType>
inline constexpr bool enable_indexed_tensor<IndexedTensor<TensorType, TensorIdxIdMapType>> = true;

template <class Type>
inline constexpr bool is_indexed_tensor_v
        = enable_indexed_tensor<std::remove_const_t<std::remove_reference_t<Type>>>;

template <class TypeSeqCharIds, class TensorType, std::size_t... I>
auto internal_index(TensorType& t, std::index_sequence<I...>)
{
    return IndexedTensor<
            TensorType,
            TensorIndexIdMap<VectorIndexIdMap<
                    ddc::type_seq_element_t<I, TypeSeqCharIds>::value,
                    typename TensorType::vector_index_set_t<I>>...>>(t);
}

template <
        class GlobalTensorIndexIdMap,
        std::size_t Is,
        class ResultIndexedTensorType,
        class... IndexedTensorType>
void internal_tensor_mul_elem(ResultIndexedTensorType& result, IndexedTensorType... t)
{
    using TensorTuple = std::tuple<typename IndexedTensorType::tensor_type...>;
    using ElementType = typename std::tuple_element_t<0, TensorTuple>::element_type;
    using Idx = get_nth_tensor_index_element_from_map_t<Is, GlobalTensorIndexIdMap>;
    if constexpr (std::is_same_v<ResultIndexedTensorType, ElementType>) {
        result
                += ((t().template get<extract_sub_tensor_element_t<
                             GlobalTensorIndexIdMap,
                             typename IndexedTensorType::index_pattern,
                             Idx>>())
                    * ...);
    } else {
        result().template get<extract_sub_tensor_element_t<
                GlobalTensorIndexIdMap,
                typename ResultIndexedTensorType::index_pattern,
                Idx>>()
                += (t().template get<extract_sub_tensor_element_t<
                            GlobalTensorIndexIdMap,
                            typename IndexedTensorType::index_pattern,
                            Idx>>()
                    * ...);
    }
}

template <
        class GlobalTensorIndexIdMap,
        class ResultIndexedTensorType,
        class... IndexedTensorType,
        std::size_t... Is>
void internal_tensor_mul(
        ResultIndexedTensorType& result,
        std::index_sequence<Is...>,
        IndexedTensorType... t)
{
    ((internal_tensor_mul_elem<GlobalTensorIndexIdMap, Is>(result, t...)), ...);
}

} // namespace tensor_tools

template <char... ids, class TensorType>
auto index(TensorType& t)
{
    using TypeSeqCharIds = ddc::detail::TypeSeq<std::integral_constant<char, ids>...>;
    return tensor_tools::internal_index<
            TypeSeqCharIds>(t, std::make_index_sequence<sizeof...(ids)>());
}

template <class... IndexedTensorType>
auto tensor_mul(IndexedTensorType... t)
{
    using GlobalTensorIndexIdMap = tensor_tools::concatenate_tensor_index_id_maps_t<
            typename IndexedTensorType::index_pattern...>;
    using TensorTuple = std::tuple<typename IndexedTensorType::tensor_type...>;
    using ElementType = typename std::tuple_element_t<0, TensorTuple>::element_type;
    static_assert(
            std::conjunction_v<std::is_same<
                    typename IndexedTensorType::tensor_type::element_type,
                    ElementType>...>,
            "All tensors must have the same element type");
    static_assert(
            std::conjunction_v<std::integral_constant<
                    bool,
                    tensor_tools::is_indexed_tensor_v<IndexedTensorType>>...>,
            "A tensor multiplication must be carried out over IndexedTensor objects");
    using ResultIndexTypeSeq = typename GlobalTensorIndexIdMap::result_indices;
    if constexpr (ddc::type_seq_size_v<ResultIndexTypeSeq> == 0) {
        ElementType result = 0.0;
        tensor_tools::internal_tensor_mul<GlobalTensorIndexIdMap>(
                result,
                std::make_index_sequence<GlobalTensorIndexIdMap::size()>(),
                t...);
        return result;
    } else {
        using ResultTensorIndexIdMap = tensor_tools::to_tensor_index_id_map_t<ResultIndexTypeSeq>;
        using ResultTensorType
                = to_tensor_t<ElementType, typename ResultTensorIndexIdMap::vector_index_sets>;
        ResultTensorType result(0.0);
        tensor_tools::IndexedTensor<ResultTensorType, ResultTensorIndexIdMap> indexed_result(
                result);
        tensor_tools::internal_tensor_mul<GlobalTensorIndexIdMap>(
                indexed_result,
                std::make_index_sequence<GlobalTensorIndexIdMap::size()>(),
                t...);
        return result;
    }
}
