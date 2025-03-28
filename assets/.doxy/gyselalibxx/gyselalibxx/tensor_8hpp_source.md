

# File tensor.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**tensor.hpp**](tensor_8hpp.md)

[Go to the documentation of this file](tensor_8hpp.md)


```C++
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

} // namespace detail

template <class ElementType, class... ValidIndexSet>
class Tensor
{
    static_assert((is_vector_index_set_v<ValidIndexSet> && ...));
    static_assert(
            (((is_covariant_vector_index_set_v<ValidIndexSet>)
              || (is_contravariant_vector_index_set_v<ValidIndexSet>))
             && ...));

public:
    using index_set = ddc::detail::TypeSeq<ValidIndexSet...>;

private:
    static constexpr std::size_t s_n_elements = (ddc::type_seq_size_v<ValidIndexSet> * ...);
    std::array<ElementType, s_n_elements> m_data;

public:
    using element_type = ElementType;

    KOKKOS_FUNCTION static constexpr std::size_t rank()
    {
        return sizeof...(ValidIndexSet);
    }

    KOKKOS_FUNCTION static constexpr std::size_t size()
    {
        return s_n_elements;
    }

private:
    static_assert(rank() > 0);

public:
    KOKKOS_DEFAULTED_FUNCTION Tensor() = default;

    explicit KOKKOS_FUNCTION Tensor(ElementType fill_value)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] = fill_value;
        }
    }

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

    template <class OElementType>
    explicit KOKKOS_FUNCTION Tensor(Tensor<OElementType, ValidIndexSet...> const& o_tensor)
        : m_data(o_tensor.m_data)
    {
    }

    template <class QueryTensorIndexElement>
    KOKKOS_FUNCTION ElementType& get()
    {
        static_assert(tensor_tools::is_tensor_index_element_v<QueryTensorIndexElement>);
        return m_data[QueryTensorIndexElement::index()];
    }

    template <class QueryTensorIndexElement>
    KOKKOS_FUNCTION ElementType const& get() const
    {
        static_assert(tensor_tools::is_tensor_index_element_v<QueryTensorIndexElement>);
        return m_data[QueryTensorIndexElement::index()];
    }

    KOKKOS_DEFAULTED_FUNCTION Tensor& operator=(Tensor const& other) = default;

    template <class... Dims>
    KOKKOS_FUNCTION Tensor& operator=(Coord<Dims...> coord)
    {
        static_assert(
                rank() == 1,
                "Copying a coordinate into a tensor is only possible for 1D tensor objects");
        static_assert(
                std::is_same_v<VectorIndexSet<Dims...>, ddc::type_seq_element_t<0, index_set>>,
                "The coordinate must have the same memory layout to make a clean conversion.");
        m_data = coord.array();
        return *this;
    }

    template <class OElementType>
    KOKKOS_FUNCTION Tensor& operator*=(OElementType val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] *= val;
        }
        return *this;
    }

    template <class OElementType>
    KOKKOS_FUNCTION Tensor& operator/=(OElementType val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] /= val;
        }
        return *this;
    }

    KOKKOS_FUNCTION Tensor& operator+=(Tensor const& val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] += val.m_data[i];
        }
        return *this;
    }

    KOKKOS_FUNCTION Tensor& operator-=(Tensor const& val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data[i] -= val.m_data[i];
        }
        return *this;
    }

    template <class OElementType>
    KOKKOS_FUNCTION Tensor operator*(OElementType val) const
    {
        Tensor result(*this);
        result *= val;
        return result;
    }

    template <class OElementType>
    KOKKOS_FUNCTION Tensor operator/(OElementType val) const
    {
        Tensor result(*this);
        result /= val;
        return result;
    }

    KOKKOS_FUNCTION Tensor operator+(Tensor const& val) const
    {
        Tensor result(*this);
        result += val;
        return result;
    }

    KOKKOS_FUNCTION Tensor operator-(Tensor const& val) const
    {
        Tensor result(*this);
        result -= val;
        return result;
    }

    KOKKOS_FUNCTION Tensor operator-() const
    {
        Tensor result;
        for (std::size_t i(0); i < s_n_elements; ++i) {
            result.m_data[i] = -m_data[i];
        }
        return result;
    }

    KOKKOS_FUNCTION bool operator==(Tensor const& o_tensor) const
    {
        bool equal(true);
        for (std::size_t i(0); i < s_n_elements; ++i) {
            equal = equal && (m_data[i] == o_tensor.m_data[i]);
        }
        return equal;
    }

    KOKKOS_FUNCTION bool operator!=(Tensor const& o_tensor) const
    {
        return !(*this == o_tensor);
    }

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

//                         Type aliases

template <class... ValidIndexSet>
using DTensor = Tensor<double, ValidIndexSet...>;

template <class ElementType, class... Dims>
using Vector = Tensor<ElementType, VectorIndexSet<Dims...>>;

template <class... Dims>
using DVector = Vector<double, Dims...>;

//                         Operators

namespace ddcHelper {

template <class... QueryIndexTag, class ElementType, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION ElementType& get(Tensor<ElementType, ValidIndexSet...>& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            ddc::detail::TypeSeq<ValidIndexSet...>,
            QueryIndexTag...>>();
}

template <class... QueryIndexTag, class ElementType, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION ElementType const& get(Tensor<ElementType, ValidIndexSet...> const& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            ddc::detail::TypeSeq<ValidIndexSet...>,
            QueryIndexTag...>>();
}

template <class ElementType, class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> to_coord(Vector<ElementType, Dims...> const& tensor)
{
    return Coord<Dims...>(get<Dims>(tensor)...);
}
} // namespace ddcHelper

template <class ElementType, class OElementType, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION Tensor<ElementType, ValidIndexSet...> operator*(
        OElementType val,
        Tensor<ElementType, ValidIndexSet...> const& tensor)
{
    return tensor * val;
}

template <class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> operator+(
        Coord<Dims...> const& coord,
        DVector<Dims...> const& tensor)
{
    return Coord<Dims...>((ddc::get<Dims>(coord) + ddcHelper::get<Dims>(tensor))...);
}

template <class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> operator-(
        Coord<Dims...> const& coord,
        DVector<Dims...> const& tensor)
{
    return Coord<Dims...>((ddc::get<Dims>(coord) - ddcHelper::get<Dims>(tensor))...);
}

template <class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...>& operator+=(
        Coord<Dims...>& coord,
        DVector<Dims...> const& tensor)
{
    ((ddc::get<Dims>(coord) += ddcHelper::get<Dims>(tensor)), ...);
    return coord;
}

template <class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...>& operator-=(
        Coord<Dims...>& coord,
        DVector<Dims...> const& tensor)
{
    ((ddc::get<Dims>(coord) -= ddcHelper::get<Dims>(tensor)), ...);
    return coord;
}
```


