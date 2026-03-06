// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "tensor_common.hpp"
#include "type_seq_tools.hpp"
#include "vector_index_tools.hpp"

/**
 * A class containing only static constexpr methods which describes
 * the identity tensor.
 * @tparam ElementType The type of the elements of the tensor (usually double/complex).
 * @tparam ValidIndexSetRow The indices of the rows.
 * @tparam ValidIndexSetCol The indices of the columns.
 */
template <class ElementType, class ValidIndexSetRow, class ValidIndexSetCol>
class IdentityTensor
{
    static_assert(is_vector_index_set_v<ValidIndexSetRow>);
    static_assert(is_vector_index_set_v<ValidIndexSetCol>);
    static_assert(ddc::type_seq_size_v<ValidIndexSetRow> == ddc::type_seq_size_v<ValidIndexSetCol>);

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
        return 2;
    }

    /**
     * @brief The size of the tensor.
     * This is the number of elements in the tensor.
     *
     * @return The number of elements in the tensor.
     */
    KOKKOS_FUNCTION static constexpr std::size_t size()
    {
        return ddc::type_seq_size_v<ValidIndexSetRow> * ddc::type_seq_size_v<ValidIndexSetCol>;
    }

    /// The TensorIndexSet describing the possible indices.
    using index_set = ddc::detail::TypeSeq<ValidIndexSetRow, ValidIndexSetCol>;

    /**
     * @brief Construct an uninitialised tensor object.
     */
    KOKKOS_DEFAULTED_FUNCTION IdentityTensor() = default;

    /**
     * @brief Get an element of the tensor.
     * @tparam QueryIndexTag A type describing the relevant index.
     * @return The relevant element of the tensor.
     */
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

    /**
     * @brief A helper type alias to get all possible indices along a
     * specified dimension.
     * @tparam Dim The dimension of interest (0 <= dim < rank()).
     */
    template <std::size_t dim>
    using vector_index_set_t = ddc::type_seq_element_t<dim, index_set>;
};

/**
 * A class containing only static constexpr methods which describes
 * the Levi-Civita tensor in Cartesian coordinates.
 * @tparam ElementType The type of the elements of the tensor (usually double/complex).
 * @tparam ValidIndexSet The indices that can be used along any of the dimensions of the tensor.
 */
template <class ElementType, class ValidIndexSet>
class CartesianLeviCivitaTensor
{
    static_assert(is_vector_index_set_v<ValidIndexSet>);
    static_assert(std::is_same_v<ValidIndexSet, vector_index_set_dual_t<ValidIndexSet>>);

private:
    static constexpr std::size_t compute_size()
   
    {
                std::size_t result = 1;
                for (std::size_t i = 0; i < rank(); ++i)
        {
                        result *= ddc::type_seq_size_v<ValidIndexSet>;
                   
        }
                return result;
           
    }

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
        return compute_size();
    }

    /// The TensorIndexSet describing the possible indices.
    using index_set = type_seq_duplicate_t<ValidIndexSet, rank()>;

    /**
     * @brief Construct an uninitialised tensor object.
     */
    KOKKOS_DEFAULTED_FUNCTION CartesianLeviCivitaTensor() = default;

    /**
     * @brief Get an element of the tensor.
     * @tparam QueryIndexTag A type describing the relevant index.
     * @return The relevant element of the tensor.
     */
    template <class QueryTensorIndexElement>
    static constexpr KOKKOS_FUNCTION ElementType get()
    {
        static_assert(tensor_tools::is_tensor_index_element_v<QueryTensorIndexElement>);
        return type_seq_permutation_parity_v<
                typename QueryTensorIndexElement::IdxTypeSeq,
                ValidIndexSet>;
    }

    /**
     * @brief A helper type alias to get all possible indices along a
     * specified dimension.
     * @tparam Dim The dimension of interest (0 <= dim < rank()).
     */
    template <std::size_t dim>
    using vector_index_set_t = ValidIndexSet;
};

/**
 * A class containing methods which describe the Levi-Civita tensor in general coordinates.
 * @tparam ElementType The type of the elements of the tensor (usually double/complex).
 * @tparam ValidIndexSet The indices that can be used along any of the dimensions of the tensor.
 */
template <class ElementType, class ValidIndexSet>
class LeviCivitaTensor
{
    static_assert(is_vector_index_set_v<ValidIndexSet>);
    static_assert(
            (is_covariant_vector_index_set_v<ValidIndexSet>)
            || (is_contravariant_vector_index_set_v<ValidIndexSet>));

private:
    ElementType m_coeff;
    static constexpr std::size_t compute_size()
    {
                std::size_t result = 1;
                for (std::size_t i = 0; i < rank(); ++i)
        {
                        result *= ddc::type_seq_size_v<ValidIndexSet>;
                   
        }
                return result;
           
    }

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
        return compute_size();
    }

    /// The TensorIndexSet describing the possible indices.
    using index_set = type_seq_duplicate_t<ValidIndexSet, rank()>;

    /**
     * @brief Construct an uninitialised tensor object.
     */
    explicit KOKKOS_FUNCTION LeviCivitaTensor(ElementType jacobian)
    {
        if constexpr (is_covariant_vector_index_set_v<ValidIndexSet>) {
            m_coeff = jacobian;
        } else {
            // Evaluating at a singular point will lead to NaN
            assert(fabs(jacobian) >= 1e-19);
            m_coeff = 1.0 / jacobian;
        }
    }

    /**
     * @brief Get an element of the tensor.
     * @tparam QueryIndexTag A type describing the relevant index.
     * @return The relevant element of the tensor.
     */
    template <class QueryTensorIndexElement>
    KOKKOS_INLINE_FUNCTION ElementType get() const
    {
        int constexpr eps = type_seq_permutation_parity_v<
                typename QueryTensorIndexElement::IdxTypeSeq,
                ValidIndexSet>;
        return eps * m_coeff;
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
 * @brief A helper function to get the value of an element of the tensor.
 * @tparam QueryIndexTag A type describing the relevant index.
 * @param tensor The tensor whose elements are examined.
 * @return The relevant element of the tensor.
 */
template <class... QueryIndexTag, class ElementType, class ValidIndexSetRow, class ValidIndexSetCol>
KOKKOS_INLINE_FUNCTION constexpr double get(
        IdentityTensor<ElementType, ValidIndexSetRow, ValidIndexSetCol> const& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            typename IdentityTensor<ElementType, ValidIndexSetRow, ValidIndexSetCol>::index_set,
            QueryIndexTag...>>();
}

/**
 * @brief A helper function to get the value of an element of the tensor.
 * @tparam QueryIndexTag A type describing the relevant index.
 * @param tensor The tensor whose elements are examined.
 * @return The relevant element of the tensor.
 */
template <class... QueryIndexTag, class ElementType, class ValidIndexSet>
KOKKOS_INLINE_FUNCTION constexpr double get(
        CartesianLeviCivitaTensor<ElementType, ValidIndexSet> const& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            typename CartesianLeviCivitaTensor<ElementType, ValidIndexSet>::index_set,
            QueryIndexTag...>>();
}

/**
 * @brief A helper function to get the value of an element of the tensor.
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

namespace detail {

template <class ElementType, class ValidIndexSetRow, class ValidIndexSetCol>
inline constexpr bool
        enable_tensor_type<IdentityTensor<ElementType, ValidIndexSetRow, ValidIndexSetCol>>
        = true;

template <class ElementType, class ValidIndexSet>
inline constexpr bool enable_tensor_type<LeviCivitaTensor<ElementType, ValidIndexSet>> = true;

template <class ElementType, class ValidIndexSet>
inline constexpr bool enable_tensor_type<CartesianLeviCivitaTensor<ElementType, ValidIndexSet>>
        = true;

} // namespace detail
