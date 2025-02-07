// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <cassert>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "tensor_index_set.hpp"
#include "vector_index_set.hpp"

namespace detail {

template <class ValidIndexSet>
struct ToCoord;

template <class... Dims>
struct ToCoord<VectorIndexSet<Dims...>>
{
    using type = Coord<Dims...>;
};

template <class ValidIndexSet>
using to_coord_t = typename ToCoord<ValidIndexSet>::type;

} // namespace detail

/**
 * @brief A class representing a Tensor.
 * @tparam ElementType The type of the elements of the tensor (usually double/complex).
 * @tparam ValidIndexSet The indices that can be used along each dimension of the tensor.
 */
template <class ElementType, class... ValidIndexSet>
class Tensor
{
    static_assert(std::conjunction_v<tensor_tools::is_vector_index_set<ValidIndexSet>...>);
    static_assert(std::conjunction_v<std::disjunction<
                          tensor_tools::is_covariant_vector_index_set<ValidIndexSet>,
                          tensor_tools::is_contravariant_vector_index_set<ValidIndexSet>>...>);

public:
    using index_set = tensor_tools::TensorIndexSet<ValidIndexSet...>;

private:
    static constexpr std::size_t s_n_elements = (ddc::type_seq_size_v<ValidIndexSet> * ...);
    std::array<ElementType, s_n_elements> m_data;

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

    using AllIndexSets = ddc::detail::TypeSeq<ValidIndexSet...>;

private:
    static_assert(rank() > 0);

public:
    /**
     * @brief Construct an uninitialised tensor object.
     */
    explicit KOKKOS_FUNCTION Tensor() {}

    /**
     * @brief Construct an uninitialised tensor object.
     */
    explicit KOKKOS_FUNCTION Tensor(ElementType fill_value)
    {
        m_data.fill(fill_value);
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
            class = std::enable_if_t<sizeof...(Params) == size()>>
    KOKKOS_FUNCTION Tensor(Params... elements) : m_data({elements...})
    {
        static_assert(
                rank() == 1,
                "Filling the tensor on initialisation is only permitted for 1D vector objects");
    }

    /**
     * @brief Construct a tensor object by copying an existing tensor.
     *
     * @param o_tensor The tensor to be copied.
     */
    template <class OElementType>
    explicit KOKKOS_FUNCTION Tensor(Tensor<OElementType, ValidIndexSet...> o_tensor)
        : m_data(o_tensor.m_data)
    {
    }

    /**
     * @brief Get a modifiable reference to an element of the tensor.
     * @tparam QueryIndexTag A type describing the relevant index.
     * @return The relevant element of the tensor.
     */
    template <class QueryTensorIndexElement>
    KOKKOS_FUNCTION ElementType& get()
    {
        static_assert(tensor_tools::is_tensor_index_element_v<QueryTensorIndexElement>);
        return m_data[QueryTensorIndexElement::index()];
    }

    /**
     * @brief Get an element of the tensor.
     * @tparam QueryIndexTag A type describing the relevant index.
     * @return The relevant element of the tensor.
     */
    template <class QueryTensorIndexElement>
    KOKKOS_FUNCTION ElementType get() const
    {
        static_assert(tensor_tools::is_tensor_index_element_v<QueryTensorIndexElement>);
        return m_data[QueryTensorIndexElement::index()];
    }

    /**
     * @brief A copy operator.
     * @param other The tensor to be copied.
     * @return A reference to the current tensor.
     */
    KOKKOS_FUNCTION Tensor& operator=(Tensor other)
    {
        for (int i(0); i < s_n_elements; ++i) {
            m_data[i] = other.m_data[i];
        }
        return *this;
    }

    /**
     * @brief An operator to multiply all the element of the current tensor by
     * a value.
     * @param val The value by which the elements should be multiplied.
     * @return A reference to the current modified tensor.
     */
    template <class OElementType>
    KOKKOS_FUNCTION Tensor& operator*=(OElementType val)
    {
        for (int i(0); i < s_n_elements; ++i) {
            m_data[i] *= val;
        }
        return *this;
    }

    /**
     * @brief An operator to add two tensors elementwise.
     * @param val The tensor that should be added to the current tensor.
     * @return A reference to the current modified tensor.
     */
    KOKKOS_FUNCTION Tensor& operator+=(Tensor val)
    {
        for (int i(0); i < s_n_elements; ++i) {
            m_data[i] += val.m_data[i];
        }
        return *this;
    }

    /**
     * @brief An operator to subtract one tensor from another elementwise.
     * @param val The tensor that should be subtracted from the current tensor.
     * @return A reference to the current modified tensor.
     */
    KOKKOS_FUNCTION Tensor& operator-=(Tensor val)
    {
        for (int i(0); i < s_n_elements; ++i) {
            m_data[i] -= val.m_data[i];
        }
        return *this;
    }

    /**
     * @brief An operator to multiply all the element of the current tensor by
     * a value.
     * @param val The value by which the elements should be multiplied.
     * @return A new tensor containing the result of the multiplication.
     */
    template <class OElementType>
    KOKKOS_FUNCTION Tensor operator*(OElementType val) const
    {
        Tensor result(*this);
        result *= val;
        return result;
    }

    /**
     * @brief An operator to add two tensors elementwise.
     * @param val The tensor that should be added to the current tensor.
     * @return A new tensor containing the result of the addition.
     */
    KOKKOS_FUNCTION Tensor operator+(Tensor val) const
    {
        Tensor result(*this);
        result += val;
        return result;
    }

    /**
     * @brief An operator to subtract one tensor from another elementwise.
     * @param val The tensor that should be subtracted from the current tensor.
     * @return A new tensor containing the result of the subtraction.
     */
    KOKKOS_FUNCTION Tensor operator-(Tensor val) const
    {
        Tensor result(*this);
        result -= val;
        return result;
    }

    /**
     * @brief A helper type alias to get all possible indices along a
     * specified dimension.
     * @tparam Dim The dimension of interest (0 <= dim < rank()).
     */
    template <std::size_t dim>
    using vector_index_set_t = ddc::type_seq_element_t<dim, AllIndexSets>;
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

//////////////////////////////////////////////////////////////////////////
//                         Operators
//////////////////////////////////////////////////////////////////////////

/**
 * @brief An operator to multiply all the element of the current tensor by
 * a value.
 * @param val The value by which the elements should be multiplied.
 * @return A new tensor containing the result of the multiplication.
 */
template <class ElementType, class OElementType, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION Tensor<ElementType, ValidIndexSet...> operator*(
        OElementType val,
        Tensor<ElementType, ValidIndexSet...> tensor)
{
    return tensor * val;
}

/**
 * An operator to add the elements of a tensor to a coordinate.
 * This can be useful in some calculations, e.g when calculating the foot
 * of a characteristic.
 * @param tensor The tensor to be converted.
 * @return The new coordinate.
 */
template <class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> operator+(Coord<Dims...> coord, DVector<Dims...> tensor)
{
    return Coord<Dims...>((ddc::get<Dims>(coord) + tensor.template get<Dims>())...);
}

/**
 * An operator to add the elements of a tensor to a coordinate.
 * This can be useful in some calculations, e.g when calculating the foot
 * of a characteristic.
 * @param tensor The tensor to be converted.
 * @return The new coordinate.
 */
template <class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> operator-(Coord<Dims...> coord, DVector<Dims...> tensor)
{
    return Coord<Dims...>((ddc::get<Dims>(coord) - tensor.template get<Dims>())...);
}

/**
 * An operator to add the elements of a tensor to a coordinate.
 * This can be useful in some calculations, e.g when calculating the foot
 * of a characteristic.
 * @param tensor The tensor to be converted.
 * @return The new coordinate.
 */
template <class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...>& operator+=(Coord<Dims...>& coord, DVector<Dims...> tensor)
{
    ((ddc::get<Dims>(coord) += tensor.template get<Dims>()), ...);
    return coord;
}

/**
 * An operator to add the elements of a tensor to a coordinate.
 * This can be useful in some calculations, e.g when calculating the foot
 * of a characteristic.
 * @param tensor The tensor to be converted.
 * @return The new coordinate.
 */
template <class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...>& operator-=(Coord<Dims...>& coord, DVector<Dims...> tensor)
{
    ((ddc::get<Dims>(coord) -= tensor.template get<Dims>()), ...);
    return coord;
}

namespace ddcHelper {

/**
 * @brief A helper function to get a modifiable reference to an element of the tensor.
 * @tparam QueryIndexTag A type describing the relevant index.
 * @param tensor The tensor whose elements are examined.
 * @return The relevant element of the tensor.
 */
template <class... QueryIndexTag, class ElementType, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION ElementType& get(Tensor<ElementType, ValidIndexSet...>& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            tensor_tools::TensorIndexSet<ValidIndexSet...>,
            QueryIndexTag...>>();
}

/**
 * @brief A helper function to get an element of the tensor.
 * @tparam QueryIndexTag A type describing the relevant index.
 * @param tensor The tensor whose elements are examined.
 * @return The relevant element of the tensor.
 */
template <class... QueryIndexTag, class ElementType, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION ElementType get(Tensor<ElementType, ValidIndexSet...> const& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            tensor_tools::TensorIndexSet<ValidIndexSet...>,
            QueryIndexTag...>>();
}

/**
 * @brief A helper function to convert a vector to a coordinate.
 * This is useful in order to add a Vector to a coordinate to obtain a
 * new coordinate (e.g. when calculating the foot of a characteristic.
 * @param tensor The tensor to be converted.
 * @return The new coordinate.
 */
template <class ElementType, class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> to_coord(Vector<ElementType, Dims...> tensor)
{
    return Coord<Dims...>(get<Dims>(tensor)...);
}
} // namespace ddcHelper
