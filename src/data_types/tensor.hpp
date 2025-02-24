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

template <class ElementType, std::size_t N>
class TupleOfElements
{
    template <typename = std::make_index_sequence<N>>
    struct impl;

    template <size_t... Is>
    struct impl<std::index_sequence<Is...>>
    {
        template <size_t>
        using wrap = ElementType;

        using type = std::tuple<wrap<Is>...>;
    };

public:
    using type = typename impl<>::type;
};

template <class TensorElementType, class OElementType>
constexpr bool is_operator_compatible_v = std::is_same_v<
        std::remove_reference_t<TensorElementType>,
        std::remove_const_t<std::remove_reference_t<OElementType>>>;

} // namespace detail

/**
 * @brief A class representing a Tensor.
 * @tparam ElementType The type of the elements of the tensor (usually double/complex).
 * @tparam ValidIndexSet The indices that can be used along each dimension of the tensor.
 */
template <class ElementType, class... ValidIndexSet>
class Tensor
{
    static_assert((is_vector_index_set_v<ValidIndexSet> && ...));
    static_assert(
            (((is_covariant_vector_index_set_v<ValidIndexSet>)
              || (is_contravariant_vector_index_set_v<ValidIndexSet>))
             && ...));

    template <class OElementType, class... OValidIndexSet>
    friend class Tensor;

public:
    /// The TensorIndexSet describing the possible indices.
    using index_set = ddc::detail::TypeSeq<ValidIndexSet...>;

private:
    static constexpr std::size_t s_n_elements = (ddc::type_seq_size_v<ValidIndexSet> * ...);
    using ArrayType = typename detail::TupleOfElements<ElementType, s_n_elements>::type;
    ArrayType m_data;

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
    template <class OElementType, size_t... Is>
    KOKKOS_FUNCTION void mul(OElementType val, std::index_sequence<Is...>)
    {
        ((std::get<Is>(m_data) *= val), ...);
    }

    template <class OElementType, size_t... Is>
    KOKKOS_FUNCTION void div(OElementType val, std::index_sequence<Is...>)
    {
        ((std::get<Is>(m_data) /= val), ...);
    }

    template <class OElementType, size_t... Is>
    KOKKOS_FUNCTION void sum(Tensor<OElementType, ValidIndexSet...> val, std::index_sequence<Is...>)
    {
        ((std::get<Is>(m_data) += std::get<Is>(val.m_data)), ...);
    }

    template <class OElementType, size_t... Is>
    KOKKOS_FUNCTION void minus(
            Tensor<OElementType, ValidIndexSet...> val,
            std::index_sequence<Is...>)
    {
        ((std::get<Is>(m_data) -= std::get<Is>(val.m_data)), ...);
    }

private:
    static_assert(rank() > 0);

public:
    /**
     * @brief Construct an uninitialised tensor object.
     */
    KOKKOS_DEFAULTED_FUNCTION Tensor() = default;

    /**
     * @brief Construct a tensor object initialised with a value.
     * @param fill_value The value with which the tensor should be filled.
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
            class = std::enable_if_t<sizeof...(Params) == size() && sizeof...(Params) != 1>>
    explicit KOKKOS_FUNCTION Tensor(Params... elements) : m_data({elements...})
    {
        static_assert(
                rank() == 1,
                "Filling the tensor on initialisation is only permitted for 1D vector objects");
    }

    /**
     * @brief Construct a 1D tensor object by providing the elements that should
     * be saved in the vector.
     *
     * @param elements The elements of the tensor.
     */
    template <
            class... Params,
            class = std::enable_if_t<
                    ((std::is_reference_v<ElementType>)&&(!std::is_const_v<ElementType>)&&(
                            std::is_convertible_v<
                                    Params,
                                    std::remove_reference_t<ElementType>> && ...))>,
            class = std::enable_if_t<sizeof...(Params) == size() && sizeof...(Params) != 1>>
    explicit KOKKOS_FUNCTION Tensor(Params&... elements) : m_data({elements...})
    {
        static_assert(std::is_reference_v<ElementType>);
        static_assert(
                rank() == 1,
                "Filling the tensor on initialisation is only permitted for 1D vector objects");
    }

    /**
     * @brief Construct a 1D tensor object from a coordinate.
     *
     * @param coord The coordinate.
     */
    template <class... Dims>
    explicit KOKKOS_FUNCTION Tensor(Coord<Dims...> coord) : m_data(coord.array())
    {
        static_assert(
                rank() == 1,
                "Filling the tensor on initialisation is only permitted for 1D vector objects");
        static_assert(
                std::is_same_v<VectorIndexSet<Dims...>, ddc::type_seq_element_t<0, index_set>>,
                "The coordinate must have the same memory layout to make a clean conversion.");
    }

    /**
     * @brief Construct a tensor object by copying an existing tensor.
     *
     * @param o_tensor The tensor to be copied.
     */
    template <class OElementType>
    explicit KOKKOS_FUNCTION Tensor(Tensor<OElementType, ValidIndexSet...> const& o_tensor)
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
        return std::get<QueryTensorIndexElement::index()>(m_data);
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
        return std::get<QueryTensorIndexElement::index()>(m_data);
    }

    /**
     * @brief A copy operator.
     * @param other The tensor to be copied.
     * @return A reference to the current tensor.
     */
    KOKKOS_DEFAULTED_FUNCTION Tensor& operator=(Tensor const& other) = default;

    /**
     * @brief A copy operator.
     * @param o_tensor The tensor to be copied.
     * @return A reference to the current tensor.
     */
    template <class OElementType>
    KOKKOS_FUNCTION Tensor& operator=(Tensor<OElementType, ValidIndexSet...> const& o_tensor)
    {
        m_data = o_tensor.m_data;
        return *this;
    }

    /**
     * @brief A copy operator.
     * @param coord The coordinate to be copied into the vector.
     * @return A reference to the current tensor.
     */
    template <class... Dims>
    KOKKOS_FUNCTION Tensor& operator=(Coord<Dims...> coord)
    {
        static_assert(
                rank() == 1,
                "Filling the tensor on initialisation is only permitted for 1D vector objects");
        static_assert(
                ddc::type_seq_same_v<
                        VectorIndexSet<Dims...>,
                        ddc::type_seq_element_t<0, index_set>>,
                "The coordinate must have the same memory layout to make a clean conversion.");
        using IndexSet = ddc::type_seq_element_t<0, index_set>;
        ((std::get<ddc::type_seq_rank_v<Dims, IndexSet>>(m_data) = ddc::get<Dims>(coord)), ...);
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
        mul(val, std::make_index_sequence<s_n_elements>());
        return *this;
    }

    /**
     * @brief An operator to divide all the element of the current tensor by
     * a value.
     * @param val The value by which the elements should be multiplied.
     * @return A reference to the current modified tensor.
     */
    template <class OElementType>
    KOKKOS_FUNCTION Tensor& operator/=(OElementType val)
    {
        div(val, std::make_index_sequence<s_n_elements>());
        return *this;
    }

    /**
     * @brief An operator to add two tensors elementwise.
     * @param val The tensor that should be added to the current tensor.
     * @return A reference to the current modified tensor.
     */
    template <
            class OElementType,
            class = std::enable_if_t<detail::is_operator_compatible_v<ElementType, OElementType>>>
    KOKKOS_FUNCTION Tensor& operator+=(Tensor<OElementType, ValidIndexSet...> const& val)
    {
        sum(val, std::make_index_sequence<s_n_elements>());
        return *this;
    }

    /**
     * @brief An operator to subtract one tensor from another elementwise.
     * @param val The tensor that should be subtracted from the current tensor.
     * @return A reference to the current modified tensor.
     */
    template <
            class OElementType,
            class = std::enable_if_t<detail::is_operator_compatible_v<ElementType, OElementType>>>
    KOKKOS_FUNCTION Tensor& operator-=(Tensor<OElementType, ValidIndexSet...> const& val)
    {
        minus(val, std::make_index_sequence<s_n_elements>());
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
     * @brief An operator to multiply all the element of the current tensor by
     * a value.
     * @param val The value by which the elements should be multiplied.
     * @return A new tensor containing the result of the multiplication.
     */
    template <class OElementType>
    KOKKOS_FUNCTION Tensor operator/(OElementType val) const
    {
        Tensor result(*this);
        result /= val;
        return result;
    }

    /**
     * @brief An operator to add two tensors elementwise.
     * @param val The tensor that should be added to the current tensor.
     * @return A new tensor containing the result of the addition.
     */
    template <
            class OElementType,
            class = std::enable_if_t<detail::is_operator_compatible_v<ElementType, OElementType>>>
    KOKKOS_FUNCTION Tensor operator+(Tensor<OElementType, ValidIndexSet...> const& val) const
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
    template <
            class OElementType,
            class = std::enable_if_t<detail::is_operator_compatible_v<ElementType, OElementType>>>
    KOKKOS_FUNCTION Tensor operator-(Tensor<OElementType, ValidIndexSet...> const& val) const
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
    using vector_index_set_t = ddc::type_seq_element_t<dim, index_set>;
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
            ddc::detail::TypeSeq<ValidIndexSet...>,
            QueryIndexTag...>>();
}

/**
 * @brief A helper function to get an element of the tensor.
 * @tparam QueryIndexTag A type describing the relevant index.
 * @param tensor The tensor whose elements are examined.
 * @return The relevant element of the tensor.
 */
template <class... QueryIndexTag, class ElementType, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION ElementType const& get(Tensor<ElementType, ValidIndexSet...> const& tensor)
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
template <class ElementType, class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> to_coord(Vector<ElementType, Dims...> const& tensor)
{
    return Coord<Dims...>(get<Dims>(tensor)...);
}
} // namespace ddcHelper

/**
 * @brief An operator to multiply all the element of the current tensor by
 * a value.
 * @param val The value by which the elements should be multiplied.
 * @return A new tensor containing the result of the multiplication.
 */
template <class ElementType, class OElementType, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION Tensor<ElementType, ValidIndexSet...> operator*(
        OElementType val,
        Tensor<ElementType, ValidIndexSet...> const& tensor)
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
KOKKOS_INLINE_FUNCTION Coord<Dims...> operator+(
        Coord<Dims...> const& coord,
        DVector<Dims...> const& tensor)
{
    return Coord<Dims...>((ddc::get<Dims>(coord) + ddcHelper::get<Dims>(tensor))...);
}

/**
 * An operator to add the elements of a tensor to a coordinate.
 * This can be useful in some calculations, e.g when calculating the foot
 * of a characteristic.
 * @param tensor The tensor to be converted.
 * @return The new coordinate.
 */
template <class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> operator-(
        Coord<Dims...> const& coord,
        DVector<Dims...> const& tensor)
{
    return Coord<Dims...>((ddc::get<Dims>(coord) - ddcHelper::get<Dims>(tensor))...);
}

/**
 * An operator to add the elements of a tensor to a coordinate.
 * This can be useful in some calculations, e.g when calculating the foot
 * of a characteristic.
 * @param tensor The tensor to be converted.
 * @return The new coordinate.
 */
template <class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...>& operator+=(
        Coord<Dims...>& coord,
        DVector<Dims...> const& tensor)
{
    ((ddc::get<Dims>(coord) += ddcHelper::get<Dims>(tensor)), ...);
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
KOKKOS_INLINE_FUNCTION Coord<Dims...>& operator-=(
        Coord<Dims...>& coord,
        DVector<Dims...> const& tensor)
{
    ((ddc::get<Dims>(coord) -= ddcHelper::get<Dims>(tensor)), ...);
    return coord;
}
