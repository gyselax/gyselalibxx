// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <cassert>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "tensor_index_tools.hpp"
#include "vector_index_tools.hpp"

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

template <class T>
inline constexpr bool enable_tensor_type = false;

} // namespace detail

template <class Type>
inline constexpr bool is_tensor_type_v
        = detail::enable_tensor_type<std::remove_const_t<std::remove_reference_t<Type>>>;

/**
 * @brief A superclass for Tensor calculations. Tensor classes containing data will inherit
 * from this class. The class Tensor will represent most Tensors but other subclasses may
 * be necessary (e.g. to access a Vector in a VectorField).
 *
 * @tparam ElementType The type of the elements of the tensor (usually double/complex).
 * @tparam LayoutType The way in which the underlying mdspan will be laid out in memory,
 *                    usually Kokkos::layout_right or Kokkos::layout_stride.
 * @tparam ValidIndexSet The indices that can be used along each dimension of the tensor.
 */
template <class ElementType, class LayoutType, class... ValidIndexSet>
class TensorCommon
{
    static_assert((is_vector_index_set_v<ValidIndexSet> && ...));
    static_assert(
            (((is_covariant_vector_index_set_v<ValidIndexSet>)
              || (is_contravariant_vector_index_set_v<ValidIndexSet>))
             && ...));

public:
    /// The TensorIndexSet describing the possible indices.
    using index_set = ddc::detail::TypeSeq<ValidIndexSet...>;

protected:
    /// The number of elements in the mdspan
    static constexpr std::size_t s_n_elements = (ddc::type_seq_size_v<ValidIndexSet> * ...);
    /// The type of the Kokkos mdspan that will be used to access the data
    using mdspan_type
            = Kokkos::mdspan<ElementType, Kokkos::extents<std::size_t, s_n_elements>, LayoutType>;
    /// The 1D object that lets us access the data
    mdspan_type m_data;

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

private:
    static_assert(rank() > 0);

protected:
    /**
     * @brief Construct an uninitialised tensor object.
     */
    KOKKOS_DEFAULTED_FUNCTION TensorCommon() = default;

    /**
     * @brief Construct a tensor object by copying an existing tensor of exactly the
     * same type. This method can be called implicitly.
     *
     * @param o_tensor The tensor to be copied.
     */
    KOKKOS_DEFAULTED_FUNCTION TensorCommon(TensorCommon const& o_tensor) = delete;

    KOKKOS_DEFAULTED_FUNCTION TensorCommon(TensorCommon&& o_tensor) = delete;

public:
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
    KOKKOS_FUNCTION ElementType const& get() const
    {
        static_assert(tensor_tools::is_tensor_index_element_v<QueryTensorIndexElement>);
        return m_data[QueryTensorIndexElement::index()];
    }

    /**
     * @brief A copy operator.
     * @param other The tensor to be copied.
     * @return A reference to the current tensor.
     */
    KOKKOS_FUNCTION TensorCommon& operator=(TensorCommon const& other)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] = other.m_data[i];
        }
        return *this;
    }

    /**
     * @brief A move assign operator.
     * @param other The tensor to be copied.
     * @return A reference to the current tensor.
     */
    KOKKOS_FUNCTION TensorCommon& operator=(TensorCommon&& other)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
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
    KOKKOS_FUNCTION TensorCommon& operator*=(OElementType val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] *= val;
        }
        return *this;
    }

    /**
     * @brief An operator to divide all the element of the current tensor by
     * a value.
     * @param val The value by which the elements should be multiplied.
     * @return A reference to the current modified tensor.
     */
    template <class OElementType>
    KOKKOS_FUNCTION TensorCommon& operator/=(OElementType val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] /= val;
        }
        return *this;
    }

    /**
     * @brief An operator to add two tensors elementwise.
     * @param val The tensor that should be added to the current tensor.
     * @return A reference to the current modified tensor.
     */
    KOKKOS_FUNCTION TensorCommon& operator+=(TensorCommon const& val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] += val.m_data[i];
        }
        return *this;
    }

    /**
     * @brief An operator to subtract one tensor from another elementwise.
     * @param val The tensor that should be subtracted from the current tensor.
     * @return A reference to the current modified tensor.
     */
    KOKKOS_FUNCTION TensorCommon& operator-=(TensorCommon const& val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] -= val.m_data[i];
        }
        return *this;
    }

    /**
     * @brief An operator to compare one tensor to another elementwise.
     * @param o_tensor The tensor that should be compared with the current tensor.
     * @return True if the tensors are equal, false otherwise.
     */
    KOKKOS_FUNCTION bool operator==(TensorCommon const& o_tensor) const
    {
        bool equal(true);
        for (std::size_t i(0); i < s_n_elements; ++i) {
            equal = equal && (m_data[i] == o_tensor.m_data[i]);
        }
        return equal;
    }

    /**
     * @brief An operator to compare one tensor to another elementwise.
     * @param o_tensor The tensor that should be compared with the current tensor.
     * @return False if the tensors are equal, true otherwise.
     */
    KOKKOS_FUNCTION bool operator!=(TensorCommon const& o_tensor) const
    {
        return !(*this == o_tensor);
    }

    /**
     * @brief A helper type alias to get all possible indices along a
     * specified dimension.
     * @tparam Dim The dimension of interest (0 <= dim < rank()).
     */
    template <std::size_t dim>
    using vector_index_set_t = ddc::type_seq_element_t<dim, index_set>;
};

//////////////////////////////////////////////////////////////////////////
//                         Operators
//////////////////////////////////////////////////////////////////////////

namespace ddcHelper {

/**
 * @brief A helper function to get a modifiable reference to an element of the tensor.
 * @tparam QueryIndexTag A type describing the relevant index.
 * @param tensor The tensor whose elements are examined.
 * @return The relevant element of the tensor.
 */
template <class... QueryIndexTag, class ElementType, class LayoutSpace, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION ElementType& get(
        TensorCommon<ElementType, LayoutSpace, ValidIndexSet...>& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            ddc::detail::TypeSeq<ValidIndexSet...>,
            QueryIndexTag...>>();
}

/**
 * @brief A helper function to get an element of the tensor.
 * @tparam QueryIndexTag A type describing the relevant index.
 * @param tensor The tensor whose elements are examined.
 * @return The relevant element of the tensor.
 */
template <class... QueryIndexTag, class ElementType, class LayoutSpace, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION ElementType const& get(
        TensorCommon<ElementType, LayoutSpace, ValidIndexSet...> const& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            ddc::detail::TypeSeq<ValidIndexSet...>,
            QueryIndexTag...>>();
}

/**
 * @brief A helper function to convert a vector to a coordinate.
 * This is useful in order to add a Vector to a coordinate to obtain a
 * new coordinate (e.g. when calculating the foot of a characteristic.
 * @param tensor The tensor to be converted.
 * @return The new coordinate.
 */
template <class ElementType, class LayoutSpace, class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> to_coord(
        TensorCommon<ElementType, LayoutSpace, ddc::detail::TypeSeq<Dims...>> const& tensor)
{
    return Coord<Dims...>(get<Dims>(tensor)...);
}
} // namespace ddcHelper

/**
 * An operator to add the elements of a tensor to a coordinate.
 * This can be useful in some calculations, e.g when calculating the foot
 * of a characteristic.
 * @param[in] coord The coordinate to which the tensor is added.
 * @param[in] tensor The tensor to be added to the coordinate.
 * @return The new coordinate.
 */
template <class... Dims, class LayoutSpace>
KOKKOS_INLINE_FUNCTION Coord<Dims...> operator+(
        Coord<Dims...> const& coord,
        TensorCommon<double, LayoutSpace, ddc::detail::TypeSeq<Dims...>> const& tensor)
{
    return Coord<Dims...>((ddc::get<Dims>(coord) + ddcHelper::get<Dims>(tensor))...);
}

/**
 * An operator to add the elements of a tensor to a coordinate.
 * This can be useful in some calculations, e.g when calculating the foot
 * of a characteristic.
 * @param[in] coord The coordinate from which the tensor is subtracted.
 * @param[in] tensor The tensor to be subtracted from the coordinate.
 * @return The new coordinate.
 */
template <class... Dims, class LayoutSpace>
KOKKOS_INLINE_FUNCTION Coord<Dims...> operator-(
        Coord<Dims...> const& coord,
        TensorCommon<double, LayoutSpace, ddc::detail::TypeSeq<Dims...>> const& tensor)
{
    return Coord<Dims...>((ddc::get<Dims>(coord) - ddcHelper::get<Dims>(tensor))...);
}

/**
 * An operator to add the elements of a tensor to a coordinate.
 * This can be useful in some calculations, e.g when calculating the foot
 * of a characteristic.
 * @param[inout] coord The coordinate to which the tensor is added.
 * @param[in] tensor The tensor to be added to the coordinate.
 * @return The new coordinate.
 */
template <class... Dims, class LayoutSpace>
KOKKOS_INLINE_FUNCTION Coord<Dims...>& operator+=(
        Coord<Dims...>& coord,
        TensorCommon<double, LayoutSpace, ddc::detail::TypeSeq<Dims...>> const& tensor)
{
    ((ddc::get<Dims>(coord) += ddcHelper::get<Dims>(tensor)), ...);
    return coord;
}

/**
 * An operator to add the elements of a tensor to a coordinate.
 * This can be useful in some calculations, e.g when calculating the foot
 * of a characteristic.
 * @param[inout] coord The coordinate from which the tensor is subtracted.
 * @param[in] tensor The tensor to be subtracted from the coordinate.
 * @return The new coordinate.
 */
template <class... Dims, class LayoutSpace>
KOKKOS_INLINE_FUNCTION Coord<Dims...>& operator-=(
        Coord<Dims...>& coord,
        TensorCommon<double, LayoutSpace, ddc::detail::TypeSeq<Dims...>> const& tensor)
{
    ((ddc::get<Dims>(coord) -= ddcHelper::get<Dims>(tensor)), ...);
    return coord;
}

/**
 * @brief An operator to multiply all a tensor by a value.
 * @param val The value by which the elements should be multiplied.
 * @param tensor The tensor being multiplied.
 * @return A new tensor containing the result of the multiplication.
 */
template <
        class TensorType,
        class OElementType,
        std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator*(OElementType val, TensorType const& tensor)
{
    TensorType result(tensor);
    result *= val;
    return result;
}

/**
 * @brief An operator to multiply all the element of the current tensor by
 * a value.
 * @param tensor The tensor being multiplied.
 * @param val The value by which the elements should be multiplied.
 * @return A new tensor containing the result of the multiplication.
 */
template <
        class TensorType,
        class OElementType,
        std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator*(TensorType const& tensor, OElementType val)
{
    TensorType result(tensor);
    result *= val;
    return result;
}

/**
 * @brief An operator to divide all the element of the current tensor by
 * a value.
 * @param tensor The tensor being divided.
 * @param val The value by which the elements should be multiplied.
 * @return A new tensor containing the result of the multiplication.
 */
template <
        class TensorType,
        class OElementType,
        std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator/(TensorType const& tensor, OElementType val)
{
    TensorType result(tensor);
    result /= val;
    return result;
}

/**
 * @brief An operator to add two tensors elementwise.
 * @param tensor The first tensor in the addition.
 * @param val The second tensor in the addition.
 * @return A new tensor containing the result of the addition.
 */
template <
        class TensorType,
        class OElementType,
        std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator+(TensorType const& tensor, OElementType val)
{
    TensorType result(tensor);
    result += val;
    return result;
}

/**
 * @brief An operator to subtract one tensor from another elementwise.
 * @param tensor The tensor which is subtracted from.
 * @param val The tensor that should be subtracted.
 * @return A new tensor containing the result of the subtraction.
 */
template <class TensorType, std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator-(TensorType const& tensor, TensorType const& val)
{
    TensorType result(tensor);
    result -= val;
    return result;
}

/**
 * @brief An operator to get the negation of a tensor elementwise.
 * @param tensor The tensor to be negated.
 * @return A new tensor containing the result of the negation.
 */
template <class TensorType, std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator-(TensorType const& tensor)
{
    TensorType result(0);
    result -= tensor;
    return result;
}
