

# File static\_tensors.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**static\_tensors.hpp**](static__tensors_8hpp.md)

[Go to the documentation of this file](static__tensors_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "type_seq_tools.hpp"
#include "vector_index_tools.hpp"

template <class ElementType, class ValidIndexSetRow, class ValidIndexSetCol>
class IdentityTensor
{
    static_assert(is_vector_index_set_v<ValidIndexSetRow>);
    static_assert(is_vector_index_set_v<ValidIndexSetCol>);
    static_assert(ddc::type_seq_size_v<ValidIndexSetRow> == ddc::type_seq_size_v<ValidIndexSetCol>);

public:
    using element_type = ElementType;

    KOKKOS_FUNCTION static constexpr std::size_t rank()
    {
        return 2;
    }

    KOKKOS_FUNCTION static constexpr std::size_t size()
    {
        return ddc::type_seq_size_v<ValidIndexSetRow> * ddc::type_seq_size_v<ValidIndexSetCol>;
    }

    using index_set = ddc::detail::TypeSeq<ValidIndexSetRow, ValidIndexSetCol>;

    KOKKOS_DEFAULTED_FUNCTION IdentityTensor() = default;

    template <class QueryTensorIndexElement>
    static constexpr KOKKOS_FUNCTION ElementType get()
    {
        static_assert(tensor_tools::is_tensor_index_element_v<QueryTensorIndexElement>);
        constexpr int row_index = ddc::type_seq_rank_v<
                typename QueryTensorIndexElement::template index_on_dim_t<0>,
                ValidIndexSetRow>;
        constexpr int col_index = ddc::type_seq_rank_v<
                typename QueryTensorIndexElement::template index_on_dim_t<1>,
                ValidIndexSetCol>;
        if constexpr (row_index == col_index) {
            return 1;
        } else {
            return 0;
        }
    }

    template <std::size_t dim>
    using vector_index_set_t = std::conditional_t<dim == 0, ValidIndexSetRow, ValidIndexSetCol>;
};

template <class ElementType, class ValidIndexSet>
class LeviCivitaTensor
{
    static_assert(is_vector_index_set_v<ValidIndexSet>);

public:
    using element_type = ElementType;

    KOKKOS_FUNCTION static constexpr std::size_t rank()
    {
        return ddc::type_seq_size_v<ValidIndexSet>;
    }

    KOKKOS_FUNCTION static constexpr std::size_t size()
    {
        return rank() * rank();
    }

    using index_set = type_seq_duplicate_t<ValidIndexSet, rank()>;

    KOKKOS_DEFAULTED_FUNCTION LeviCivitaTensor() = default;

    template <class QueryTensorIndexElement>
    static constexpr KOKKOS_FUNCTION ElementType get()
    {
        static_assert(tensor_tools::is_tensor_index_element_v<QueryTensorIndexElement>);
        if constexpr (type_seq_has_unique_elements_v<
                              typename QueryTensorIndexElement::IdxTypeSeq>) {
            return type_seq_permutation_parity_v<
                    typename QueryTensorIndexElement::IdxTypeSeq,
                    ValidIndexSet>;
        } else {
            return 0;
        }
    }

    template <std::size_t dim>
    using vector_index_set_t = ValidIndexSet;
};


namespace ddcHelper {
template <class... QueryIndexTag, class ElementType, class ValidIndexSetRow, class ValidIndexSetCol>
KOKKOS_INLINE_FUNCTION constexpr double get(
        IdentityTensor<ElementType, ValidIndexSetRow, ValidIndexSetCol> const& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            typename IdentityTensor<ElementType, ValidIndexSetRow, ValidIndexSetCol>::index_set,
            QueryIndexTag...>>();
}

template <class... QueryIndexTag, class ElementType, class ValidIndexSet>
KOKKOS_INLINE_FUNCTION constexpr double get(
        LeviCivitaTensor<ElementType, ValidIndexSet> const& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            typename LeviCivitaTensor<ElementType, ValidIndexSet>::index_set,
            QueryIndexTag...>>();
}
} // namespace ddcHelper
```


