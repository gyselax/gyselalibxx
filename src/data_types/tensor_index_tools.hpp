// SPDX-License-Identifier: MIT
#pragma once
#include <tuple>

#include "tensor_index_tools.hpp"
#include "vector_index_set.hpp"

namespace tensor_tools {

template <class... ValidIndexSet>
class TensorIndexSet;

template <class ValidatingTensorIndexSet, class... Dims>
class TensorIndexElement;

template <class... ValidIndexSet>
class TensorIndexSet
{
    static_assert(std::conjunction_v<is_vector_index_set<ValidIndexSet>...>);

private:
    static constexpr std::size_t s_n_elements = (ddc::type_seq_size_v<ValidIndexSet> * ...);

    using AllIndexSets = ddc::detail::TypeSeq<ValidIndexSet...>;

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
        return sizeof...(ValidIndexSet);
    }

    /**
     * @brief The size of the tensor.
     * This is the number of elements in the tensor.
     *
     * @return The number of elements in the tensor.
     */
    KOKKOS_FUNCTION static constexpr std::size_t size()
    {
        return s_n_elements;
    }

    template <std::size_t Idx>
    static constexpr std::size_t get_stride()
    {
        static_assert(Idx < rank());
        if constexpr (Idx >= rank() - 1) {
            return 1;
        } else {
            std::size_t shape_elem
                    = ddc::type_seq_size_v<ddc::type_seq_element_t<Idx + 1, AllIndexSets>>;
            return shape_elem * get_stride<Idx + 1>();
        }
    }

    template <std::size_t IDim>
    using get_vector_index_set_along_dim_t = ddc::type_seq_element_t<IDim, AllIndexSets>;

private:
    template <class TensorIndexSetType, class... Dims>
    friend class TensorIndexElement;
};

template <class ValidatingTensorIndexSet, class... Dims>
class TensorIndexElement
{
    static_assert(is_tensor_index_set_v<ValidatingTensorIndexSet>);
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
    static_assert(rank() == ValidatingTensorIndexSet::rank(), "Wrong number of indices provided");
    template <std::size_t... Is>
    static constexpr bool valid_indices(std::index_sequence<Is...>)
    {
        return ((ddc::in_tags_v<
                 ddc::type_seq_element_t<Is, IdxTypeSeq>,
                 typename ValidatingTensorIndexSet::get_vector_index_set_along_dim_t<Is>>)&&...);
    }

    static_assert(valid_indices(std::make_index_sequence<sizeof...(Dims)>()));

    template <std::size_t... Is>
    KOKKOS_FUNCTION static constexpr std::size_t internal_index(std::index_sequence<Is...>)
    {
        return ((ddc::type_seq_rank_v<
                         ddc::type_seq_element_t<Is, IdxTypeSeq>,
                         typename ValidatingTensorIndexSet::get_vector_index_set_along_dim_t<
                                 Is>> * ValidatingTensorIndexSet::template get_stride<Is>())
                + ...);
    }

public:
    KOKKOS_FUNCTION static constexpr std::size_t index()
    {
        return internal_index(std::make_index_sequence<rank()>());
    }

    template <std::size_t dim_idx>
    using index_on_dim_t = ddc::type_seq_element_t<dim_idx, IdxTypeSeq>;

    using tensor_index_set = ValidatingTensorIndexSet;
};

namespace details {
template <class T>
inline constexpr bool enable_tensor_index_set = false;

template <class... ValidIndexSet>
inline constexpr bool enable_tensor_index_set<TensorIndexSet<ValidIndexSet...>> = true;

template <class T>
inline constexpr bool enable_tensor_index_element = false;

template <class ValidatingTensorIndexSet, class... Dims>
inline constexpr bool
        enable_tensor_index_element<TensorIndexElement<ValidatingTensorIndexSet, Dims...>> = true;

template <class TypeSeqValidIndexSet>
struct ToTensorIdxSet;

template <class... ValidIndexSet>
struct ToTensorIdxSet<ddc::detail::TypeSeq<ValidIndexSet...>>
{
    using type = TensorIndexSet<typename ValidIndexSet::possible_idx_values...>;
};

template <class ValidatingTensorIndexSet, class TypeSeqValidIndexSet>
struct ToTensorIdxElement;

template <class ValidatingTensorIndexSet, class... ValidIndexSet>
struct ToTensorIdxElement<ValidatingTensorIndexSet, ddc::detail::TypeSeq<ValidIndexSet...>>
{
    using type = TensorIndexElement<ValidatingTensorIndexSet, ValidIndexSet...>;
};
}

template <typename Type>
inline constexpr bool is_tensor_index_set_v
        = details::enable_tensor_index_set<std::remove_const_t<std::remove_reference_t<Type>>>;

template <class TypeSeqValidIndexSet>
using to_tensor_idx_set_t = typename ToTensorIdxSet<TypeSeqValidIndexSet>::type;

template <typename Type>
inline constexpr bool is_tensor_index_element_v
        = details::enable_tensor_index_element<std::remove_const_t<std::remove_reference_t<Type>>>;

template <class ValidatingTensorIndexSet, class TypeSeqValidIndexSet>
using to_tensor_idx_element_t =
        typename ToTensorIdxElement<ValidatingTensorIndexSet, TypeSeqValidIndexSet>::type;

