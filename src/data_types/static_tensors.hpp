// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "type_seq_tools.hpp"
#include "vector_index_tools.hpp"

/**
 * A class containing only static constexpr methods which describes
 * the Levi-Civita tensor.
 * @tparam ElementType The type of the elements of the tensor (usually double/complex).
 * @tparam ValidIndexSet The indices that can be used along any of the dimensions of the tensor.
 */
template <class ElementType, class ValidIndexSet>
class LeviCivitaTensor
{
    static_assert(is_vector_index_set_v<ValidIndexSet>);
    static_assert(is_covariant_vector_index_set_v<ValidIndexSet>,
            "The Levi-Civita tensor is defined from covariant components");

public:
    /// The type of the elements of the tensor.
    using element_type = ElementType;

    /**
     * @brief The rank of the tensor.
     * This is equivalent to the number of indices required to
     * access an element of the tensor.
     *
     * @return The rank of the tensor.
     */
    KOKKOS_FUNCTION static constexpr std::size_t rank()
    {
        return ddc::type_seq_size_v<ValidIndexSet>;
    }

    /**
     * @brief The size of the tensor.
     * This is the number of elements in the tensor.
     *
     * @return The number of elements in the tensor.
     */
    KOKKOS_FUNCTION static constexpr std::size_t size()
    {
        return rank() * rank();
    }

    /// The TensorIndexSet describing the possible indices.
    using index_set = type_seq_duplicate_t<ValidIndexSet, rank()>;

    /**
     * @brief Construct an uninitialised tensor object.
     */
    KOKKOS_DEFAULTED_FUNCTION LeviCivitaTensor() = default;

    /**
     * @brief Get an element of the tensor.
     * @tparam QueryIndexTag A type describing the relevant index.
     * @return The relevant element of the tensor.
     */
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

    /**
     * @brief A helper type alias to get all possible indices along a
     * specified dimension.
     * @tparam Dim The dimension of interest (0 <= dim < rank()).
     */
    template <std::size_t dim>
    using vector_index_set_t = ValidIndexSet;
};


namespace ddcHelper {
/**
 * @brief A helper function to get a modifiable reference to an element of the tensor.
 * @tparam QueryIndexTag A type describing the relevant index.
 * @param tensor The tensor whose elements are examined.
 * @return The relevant element of the tensor.
 */
template <class... QueryIndexTag, class ElementType, class ValidIndexSet>
KOKKOS_INLINE_FUNCTION constexpr double get(
        LeviCivitaTensor<ElementType, ValidIndexSet> const& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            typename LeviCivitaTensor<ElementType, ValidIndexSet>::index_set,
            QueryIndexTag...>>();
}
} // namespace ddcHelper
