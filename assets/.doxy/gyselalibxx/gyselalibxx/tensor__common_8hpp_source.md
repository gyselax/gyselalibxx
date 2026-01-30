

# File tensor\_common.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**tensor\_common.hpp**](tensor__common_8hpp.md)

[Go to the documentation of this file](tensor__common_8hpp.md)


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

template <class T>
inline constexpr bool enable_tensor_type = false;

} // namespace detail

template <class Type>
inline constexpr bool is_tensor_type_v
        = detail::enable_tensor_type<std::remove_const_t<std::remove_reference_t<Type>>>;

template <class DataStorageType, class... ValidIndexSet>
class TensorCommon
{
    static_assert((is_vector_index_set_v<ValidIndexSet> && ...));
    static_assert(
            (((is_covariant_vector_index_set_v<ValidIndexSet>)
              || (is_contravariant_vector_index_set_v<ValidIndexSet>))
             && ...));

public:
    using index_set = ddc::detail::TypeSeq<ValidIndexSet...>;

protected:
    static constexpr std::size_t s_n_elements = (ddc::type_seq_size_v<ValidIndexSet> * ...);

    DataStorageType m_data;

public:
    using element_type = typename DataStorageType::element_type;

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

protected:
    KOKKOS_DEFAULTED_FUNCTION TensorCommon() = default;

    KOKKOS_DEFAULTED_FUNCTION TensorCommon(TensorCommon const& o_tensor) = default;

    KOKKOS_DEFAULTED_FUNCTION TensorCommon(TensorCommon&& o_tensor) = default;

public:
    template <class QueryTensorIndexElement>
    KOKKOS_FUNCTION element_type& get()
    {
        static_assert(tensor_tools::is_tensor_index_element_v<QueryTensorIndexElement>);
        typename DataStorageType::mdspan_type data = m_data();
        return data[QueryTensorIndexElement::index()];
    }

    template <class QueryTensorIndexElement>
    KOKKOS_FUNCTION element_type const& get() const
    {
        static_assert(tensor_tools::is_tensor_index_element_v<QueryTensorIndexElement>);
        typename DataStorageType::const_mdspan_type data = m_data();
        return data[QueryTensorIndexElement::index()];
    }

    KOKKOS_FUNCTION TensorCommon& operator=(TensorCommon const& other)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data()[i] = other.m_data()[i];
        }
        return *this;
    }

    KOKKOS_FUNCTION TensorCommon& operator=(TensorCommon&& other)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data()[i] = other.m_data()[i];
        }
        return *this;
    }

    template <class Oelement_type>
    KOKKOS_FUNCTION TensorCommon& operator*=(Oelement_type val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data()[i] *= val;
        }
        return *this;
    }

    template <class Oelement_type>
    KOKKOS_FUNCTION TensorCommon& operator/=(Oelement_type val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data()[i] /= val;
        }
        return *this;
    }

    KOKKOS_FUNCTION TensorCommon& operator+=(TensorCommon const& val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data()[i] += val.m_data()[i];
        }
        return *this;
    }

    KOKKOS_FUNCTION TensorCommon& operator-=(TensorCommon const& val)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data()[i] -= val.m_data()[i];
        }
        return *this;
    }

    KOKKOS_FUNCTION bool operator==(TensorCommon const& o_tensor) const
    {
        bool equal(true);
        for (std::size_t i(0); i < s_n_elements; ++i) {
            equal = equal && (m_data()[i] == o_tensor.m_data()[i]);
        }
        return equal;
    }

    KOKKOS_FUNCTION bool operator!=(TensorCommon const& o_tensor) const
    {
        return !(*this == o_tensor);
    }

    template <std::size_t dim>
    using vector_index_set_t = ddc::type_seq_element_t<dim, index_set>;
};

//                         Operators

namespace detail {
template <
        class InputStorageType,
        class OutputStorageType,
        class... OutputValidIndexSet,
        class... InputValidIndexSet,
        std::size_t... InternalIdx>
KOKKOS_INLINE_FUNCTION void assign_elements(
        TensorCommon<OutputStorageType, OutputValidIndexSet...>& tensor_to_fill,
        TensorCommon<InputStorageType, InputValidIndexSet...> const& tensor_input,
        std::index_sequence<InternalIdx...>)
{
    using InputTensorIndexSet = ddc::detail::TypeSeq<InputValidIndexSet...>;
    ((tensor_to_fill.template get<tensor_tools::to_tensor_index_element_t<
              ddc::detail::TypeSeq<OutputValidIndexSet...>,
              typename tensor_tools::get_nth_tensor_index_element_t<
                      InternalIdx,
                      InputTensorIndexSet>::IdxTypeSeq>>()
      = tensor_input.template get<
              tensor_tools::get_nth_tensor_index_element_t<InternalIdx, InputTensorIndexSet>>()),
     ...);
}
} // namespace detail

namespace ddcHelper {

template <class... QueryIndexTag, class storage_type, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION typename storage_type::element_type& get(
        TensorCommon<storage_type, ValidIndexSet...>& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            ddc::detail::TypeSeq<ValidIndexSet...>,
            QueryIndexTag...>>();
}

template <class... QueryIndexTag, class storage_type, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION typename storage_type::element_type const& get(
        TensorCommon<storage_type, ValidIndexSet...> const& tensor)
{
    return tensor.template get<tensor_tools::TensorIndexElement<
            ddc::detail::TypeSeq<ValidIndexSet...>,
            QueryIndexTag...>>();
}

template <class storage_type, class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> to_coord(
        TensorCommon<storage_type, ddc::detail::TypeSeq<Dims...>> const& tensor)
{
    return Coord<Dims...>(get<Dims>(tensor)...);
}

template <
        class InputStorageType,
        class OutputStorageType,
        class... OutputValidIndexSet,
        class... InputValidIndexSet,
        std::size_t... InternalIdx>
KOKKOS_INLINE_FUNCTION void assign_elements(
        TensorCommon<OutputStorageType, OutputValidIndexSet...>& tensor_to_fill,
        TensorCommon<InputStorageType, InputValidIndexSet...> const& tensor_input)
{
    detail::assign_elements(
            tensor_to_fill,
            tensor_input,
            std::make_index_sequence<
                    TensorCommon<InputStorageType, InputValidIndexSet...>::size()>());
}

} // namespace ddcHelper

template <class storage_type, class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> operator+(
        Coord<Dims...> const& coord,
        TensorCommon<storage_type, ddc::detail::TypeSeq<Dims...>> const& tensor)
{
    static_assert(std::is_same_v<typename storage_type::element_type, double>);
    return Coord<Dims...>((ddc::get<Dims>(coord) + ddcHelper::get<Dims>(tensor))...);
}

template <class storage_type, class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...> operator-(
        Coord<Dims...> const& coord,
        TensorCommon<storage_type, ddc::detail::TypeSeq<Dims...>> const& tensor)
{
    static_assert(std::is_same_v<typename storage_type::element_type, double>);
    return Coord<Dims...>((ddc::get<Dims>(coord) - ddcHelper::get<Dims>(tensor))...);
}

template <class storage_type, class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...>& operator+=(
        Coord<Dims...>& coord,
        TensorCommon<storage_type, ddc::detail::TypeSeq<Dims...>> const& tensor)
{
    static_assert(std::is_same_v<typename storage_type::element_type, double>);
    ((ddc::get<Dims>(coord) += ddcHelper::get<Dims>(tensor)), ...);
    return coord;
}

template <class storage_type, class... Dims>
KOKKOS_INLINE_FUNCTION Coord<Dims...>& operator-=(
        Coord<Dims...>& coord,
        TensorCommon<storage_type, ddc::detail::TypeSeq<Dims...>> const& tensor)
{
    static_assert(std::is_same_v<typename storage_type::element_type, double>);
    ((ddc::get<Dims>(coord) -= ddcHelper::get<Dims>(tensor)), ...);
    return coord;
}

template <
        class TensorType,
        class Oelement_type,
        std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator*(Oelement_type val, TensorType const& tensor)
{
    TensorType result(tensor);
    result *= val;
    return result;
}

template <
        class TensorType,
        class Oelement_type,
        std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator*(TensorType const& tensor, Oelement_type val)
{
    TensorType result(tensor);
    result *= val;
    return result;
}

template <
        class TensorType,
        class Oelement_type,
        std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator/(TensorType const& tensor, Oelement_type val)
{
    TensorType result(tensor);
    result /= val;
    return result;
}

template <
        class TensorType,
        class Oelement_type,
        std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator+(TensorType const& tensor, Oelement_type val)
{
    TensorType result(tensor);
    result += val;
    return result;
}

template <class TensorType, std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator-(TensorType const& tensor, TensorType const& val)
{
    TensorType result(tensor);
    result -= val;
    return result;
}

template <class TensorType, std::enable_if_t<is_tensor_type_v<TensorType>, bool> = true>
KOKKOS_FUNCTION TensorType operator-(TensorType const& tensor)
{
    TensorType result(0);
    result -= tensor;
    return result;
}
```


