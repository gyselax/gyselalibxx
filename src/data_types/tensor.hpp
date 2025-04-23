// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <cassert>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "tensor_common.hpp"
#include "tensor_index_tools.hpp"
#include "vector_index_tools.hpp"

/**
 * @brief A class representing a Tensor.
 * @tparam ElementType The type of the elements of the tensor (usually double/complex).
 * @tparam ValidIndexSet The indices that can be used along each dimension of the tensor.
 */
template <class ElementType, class... ValidIndexSet>
class Tensor : public TensorCommon<ElementType, detail::LocalExecutionSpace, ValidIndexSet...>
{
    using base_type = TensorCommon<ElementType, detail::LocalExecutionSpace, ValidIndexSet...>;

public:
    /// The TensorIndexSet describing the possible indices.
    using typename base_type::index_set;

    using base_type::rank;

private:
    using base_type::s_n_elements;
    using typename base_type::mdspan_type;
    std::array<ElementType, s_n_elements> m_data_alloc;
    using base_type::m_data;

public:
    /**
     * @brief Construct an uninitialised tensor object.
     */
    KOKKOS_FUNCTION Tensor() : base_type(mdspan_type(m_data_alloc.data())) {}

    /**
     * @brief Construct a tensor object initialised with a value.
     * @param fill_value The value with which the tensor should be filled.
     */
    explicit KOKKOS_FUNCTION Tensor(ElementType fill_value)
        : base_type(mdspan_type(m_data_alloc.data()))
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] = fill_value;
        }
    }

    /**
     * @brief Construct a 1D tensor object by providing the elements that should
     * be saved in the vector.
     *
     * @param elements The elements of the tensor.
     */
    template <
            class... Params,
            class = std::enable_if_t<(std::is_convertible_v<Params, ElementType> && ...)>,
            class = std::enable_if_t<
                    sizeof...(Params) == base_type::size() && sizeof...(Params) != 1>>
    explicit KOKKOS_FUNCTION Tensor(Params... elements)
        : base_type(mdspan_type(m_data_alloc.data()))
    {
        static_assert(
                rank() == 1,
                "Filling the tensor on initialisation is only permitted for 1D vector objects");
        m_data_alloc = std::array<ElementType, base_type::size()>({elements...});
    }

    /**
     * @brief Construct a 1D tensor object from a coordinate.
     *
     * @param coord The coordinate.
     */
    template <class... Dims>
    explicit KOKKOS_FUNCTION Tensor(Coord<Dims...> coord)
        : base_type(mdspan_type(m_data_alloc.data()))
    {
        static_assert(
                rank() == 1,
                "Filling the tensor on initialisation is only permitted for 1D vector objects");
        static_assert(
                std::is_same_v<VectorIndexSet<Dims...>, ddc::type_seq_element_t<0, index_set>>,
                "The coordinate must have the same memory layout to make a clean conversion.");
        m_data_alloc = coord.array();
    }

    /**
     * @brief Construct a tensor object by copying an existing compatible tensor.
     * A tensor is compatible if it is defined on the same dimensions.
     *
     * @param o_tensor The tensor to be copied.
     */
    template <class OTensorType, std::enable_if_t<is_tensor_type_v<OTensorType>, bool> = true>
    explicit KOKKOS_FUNCTION Tensor(const OTensorType& o_tensor)
        : base_type(mdspan_type(m_data_alloc.data()))
    {
        static_assert(
                std::is_same_v<typename OTensorType::index_set, index_set>,
                "The coordinate must have the same memory layout to make a clean conversion.");
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] = o_tensor.m_data[i];
        }
    }

    /**
     * @brief Construct a tensor object by copying an existing tensor of exactly the
     * same type. This method can be called implicitly.
     *
     * @param o_tensor The tensor to be copied.
     */
    KOKKOS_FUNCTION Tensor(Tensor const& o_tensor) : base_type(mdspan_type(m_data_alloc.data()))
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] = o_tensor.m_data[i];
        }
    }

    /**
     * @brief Construct a tensor object by moving an existing tensor of exactly the
     * same type. This method can be called implicitly.
     *
     * @param o_tensor The tensor to be copied.
     */
    explicit KOKKOS_FUNCTION Tensor(Tensor&& o_tensor)
        : m_data_alloc(std::move(o_tensor.m_data_alloc))
        , base_type(mdspan_type(m_data_alloc.data()))
    {
    }

    /**
     * @brief A copy operator.
     * @param other The tensor to be copied.
     * @return A reference to the current tensor.
     */
    KOKKOS_FUNCTION Tensor& operator=(Tensor const& other)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] = other.m_data[i];
        }
        return *this;
    }
};

namespace detail {

template <class ElementType, class TypeSeqValidIndexSet>
struct ToTensor;

template <class ElementType, class... ValidIndexSet>
struct ToTensor<ElementType, ddc::detail::TypeSeq<ValidIndexSet...>>
{
    using type = Tensor<ElementType, ValidIndexSet...>;
};

} // namespace detail

template <class ElementType, class TypeSeqValidIndexSet>
using to_tensor_t = typename detail::ToTensor<ElementType, TypeSeqValidIndexSet>::type;

//////////////////////////////////////////////////////////////////////////
//                         Type aliases
//////////////////////////////////////////////////////////////////////////

/**
 * @brief A helper type alias to get a tensor containing doubles.
 * @tparam ValidIndexSet The indices that can be used along each dimension of the tensor.
 */
template <class... ValidIndexSet>
using DTensor = Tensor<double, ValidIndexSet...>;

/**
 * @brief A helper type alias to get a 1D tensor (a vector).
 * @tparam ElementType The type of the elements of the tensor (usually double/complex).
 * @tparam Dims The dimensions that can be used to index the vector.
 */
template <class ElementType, class... Dims>
using Vector = Tensor<ElementType, VectorIndexSet<Dims...>>;

/**
 * @brief A helper type alias to get a 1D tensor (a vector) of doubles.
 * @tparam Dims The dimensions that can be used to index the vector.
 */
template <class... Dims>
using DVector = Vector<double, Dims...>;

namespace detail {

template <class ElementType, class... ValidIndexSet>
inline constexpr bool enable_tensor_type<Tensor<ElementType, ValidIndexSet...>> = true;

}
