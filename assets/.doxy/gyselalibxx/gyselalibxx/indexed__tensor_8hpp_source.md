

# File indexed\_tensor.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**indexed\_tensor.hpp**](indexed__tensor_8hpp.md)

[Go to the documentation of this file](indexed__tensor_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "tensor.hpp"
#include "tensor_index_tools.hpp"
#include "vector_index_tools.hpp"

namespace tensor_tools {

template <class T>
inline constexpr bool enable_indexed_tensor = false;

template <class TensorType, class TypeSeqVectorIndexIdMap>
class IndexedTensor
{
public:
    using index_pattern = TypeSeqVectorIndexIdMap;
    using tensor_type = TensorType;

    static_assert(
            type_seq_has_unique_elements_v<index_pattern>,
            "You should not have more than two of any one index in an index expression. "
            "Additionally repeated indices should not be associated with two covariant or two "
            "contravariant indices.");

private:
    TensorType& m_tensor;

public:
    explicit KOKKOS_FUNCTION IndexedTensor(TensorType& tensor) : m_tensor(tensor) {}

    IndexedTensor(IndexedTensor const&) = delete;

    IndexedTensor(IndexedTensor&&) = delete;

    KOKKOS_DEFAULTED_FUNCTION ~IndexedTensor() noexcept = default;

    IndexedTensor& operator=(IndexedTensor const&) = delete;

    IndexedTensor& operator=(IndexedTensor&&) = delete;

    KOKKOS_FUNCTION TensorType& operator()()
    {
        return m_tensor;
    }

    KOKKOS_FUNCTION TensorType const& operator()() const
    {
        return m_tensor;
    }
};

template <class TensorType, class TypeSeqVectorIndexIdMap>
inline constexpr bool
        enable_indexed_tensor<IndexedTensor<TensorType, TypeSeqVectorIndexIdMap>> = true;

template <class Type>
inline constexpr bool is_indexed_tensor_v
        = enable_indexed_tensor<std::remove_const_t<std::remove_reference_t<Type>>>;

namespace details {
template <class TypeSeqCharIds, class TensorType, std::size_t... I>
KOKKOS_FUNCTION auto internal_index(TensorType const& tensor, std::index_sequence<I...>)
{
    return IndexedTensor<
            const TensorType,
            ddc::detail::TypeSeq<VectorIndexIdMap<
                    ddc::type_seq_element_t<I, TypeSeqCharIds>::value,
                    typename TensorType::template vector_index_set_t<I>>...>>(tensor);
}

template <
        class GlobalTensorIndexIdMap,
        std::size_t Is,
        class ResultIndexedTensorType,
        class... IndexedTensorType>
KOKKOS_FUNCTION void internal_tensor_mul_elem(
        ResultIndexedTensorType& result,
        IndexedTensorType const&... t)
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
KOKKOS_FUNCTION void internal_tensor_mul(
        ResultIndexedTensorType& result,
        std::index_sequence<Is...>,
        IndexedTensorType const&... t)
{
    ((internal_tensor_mul_elem<GlobalTensorIndexIdMap, Is>(result, t...)), ...);
}

template <char ID, class AllIndexIdMaps>
KOKKOS_INLINE_FUNCTION void check_id_validity()
{
    constexpr std::size_t n_occurrences
            = char_occurrences_v<ID, index_identifiers_t<AllIndexIdMaps>>;
    if constexpr (n_occurrences == 2) {
        using RelevantVectorIndexSets = relevant_vector_index_sets_t<ID, AllIndexIdMaps>;
        using Set1 = ddc::type_seq_element_t<0, RelevantVectorIndexSets>;
        using Set2 = ddc::type_seq_element_t<1, RelevantVectorIndexSets>;
        constexpr bool has_covariant_idx = (is_covariant_vector_index_set_v<Set1>)
                                           || (is_covariant_vector_index_set_v<Set2>);
        constexpr bool has_contravariant_idx = (is_contravariant_vector_index_set_v<Set1>)
                                               || (is_contravariant_vector_index_set_v<Set2>);
        static_assert(
                has_covariant_idx && has_contravariant_idx,
                "Repeated indices should not be associated with two covariant or two contravariant "
                "indices.");
        static_assert(
                std::is_same_v<get_contravariant_dims_t<Set1>, get_contravariant_dims_t<Set2>>,
                "Cannot sum over incompatible VectorIndexSets.");
    } else {
        static_assert(
                n_occurrences == 1,
                "You should not have more than two of any one index in an index expression.");
    }
}

template <class AllIndexIdMaps, class UniqueIds, std::size_t... Is>
KOKKOS_INLINE_FUNCTION void check_all_id_validity(std::index_sequence<Is...>)
{
    ((check_id_validity<ddc::type_seq_element_t<Is, UniqueIds>::value, AllIndexIdMaps>()), ...);
}

} // namespace details

} // namespace tensor_tools

template <char... ids, class TensorType>
KOKKOS_FUNCTION auto index(TensorType const& tensor)
{
    static_assert(
            sizeof...(ids) == TensorType::rank(),
            "The number of indices provided must be equal to the rank of the tensor.");
    using namespace tensor_tools;
    using TypeSeqCharIds = ddc::detail::TypeSeq<std::integral_constant<char, ids>...>;
    return details::internal_index<
            TypeSeqCharIds>(tensor, std::make_index_sequence<sizeof...(ids)>());
}

template <class... IndexedTensorType>
KOKKOS_FUNCTION auto tensor_mul(IndexedTensorType... tensor_to_mul)
{
    using namespace tensor_tools;
    static_assert(
            (is_indexed_tensor_v<IndexedTensorType> && ...),
            "A tensor multiplication must be carried out over IndexedTensor objects");
    // Get the TypeSeq of VectorIndexIdMaps describing all the indices which appear in the calculation.
    using AllIndexIdMaps = type_seq_cat_t<typename IndexedTensorType::index_pattern...>;
    using UniqueIndexIds = type_seq_unique_t<index_identifiers_t<AllIndexIdMaps>>;
    details::check_all_id_validity<AllIndexIdMaps, UniqueIndexIds>(
            std::make_index_sequence<ddc::type_seq_size_v<UniqueIndexIds>>());
    // Use the first tensor argument to extract the element type of the result.
    using ElementType = std::tuple_element_t<
            0,
            std::tuple<typename IndexedTensorType::tensor_type::element_type...>>;
    static_assert(
            (std::is_same_v<
                     typename IndexedTensorType::tensor_type::element_type,
                     ElementType> && ...),
            "All tensors in a tensor_mul must have the same element type.");
    using ResultIndexTypeSeq = non_repeated_indices_t<AllIndexIdMaps>;
    if constexpr (ddc::type_seq_size_v<ResultIndexTypeSeq> == 0) {
        ElementType result = 0.0;
        details::internal_tensor_mul<AllIndexIdMaps>(
                result,
                std::make_index_sequence<
                        details::CalculateSize<unique_indices_t<AllIndexIdMaps>>::value>(),
                tensor_to_mul...);
        return result;
    } else {
        using ResultTensorType
                = to_tensor_t<ElementType, get_type_seq_vector_index_set_t<ResultIndexTypeSeq>>;
        ResultTensorType result(0.0);
        IndexedTensor<ResultTensorType, ResultIndexTypeSeq> indexed_result(result);
        details::internal_tensor_mul<AllIndexIdMaps>(
                indexed_result,
                std::make_index_sequence<
                        details::CalculateSize<unique_indices_t<AllIndexIdMaps>>::value>(),
                tensor_to_mul...);
        return result;
    }
}
```


