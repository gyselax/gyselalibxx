

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
#include "tensor_common.hpp"
#include "tensor_index_tools.hpp"
#include "vector_index_tools.hpp"

namespace detail {
template <class ElementType, std::size_t n_elements>
struct TensorDataInnards
{
    using element_type = ElementType;

    static constexpr std::size_t s_n_elements = n_elements;

    using mdspan_type = Kokkos::
            mdspan<ElementType, Kokkos::extents<std::size_t, s_n_elements>, Kokkos::layout_right>;
    using const_mdspan_type = Kokkos::mdspan<
            const ElementType,
            Kokkos::extents<std::size_t, s_n_elements>,
            Kokkos::layout_right>;

    std::array<ElementType, s_n_elements> m_data_alloc;

    KOKKOS_INLINE_FUNCTION const_mdspan_type operator()() const
    {
        return const_mdspan_type(m_data_alloc.data());
    }

    KOKKOS_INLINE_FUNCTION mdspan_type operator()()
    {
        return mdspan_type(m_data_alloc.data());
    }
};
} // namespace detail


template <class ElementType, class ValidIndexSetFirstDim, class... ValidIndexSet>
class Tensor
    : public TensorCommon<
              detail::TensorDataInnards<
                      ElementType,
                      (ddc::type_seq_size_v<ValidIndexSet> * ...
                       * ddc::type_seq_size_v<ValidIndexSetFirstDim>)>,
              ValidIndexSetFirstDim,
              ValidIndexSet...>
{
    using base_type = TensorCommon<
            detail::TensorDataInnards<
                    ElementType,
                    (ddc::type_seq_size_v<ValidIndexSet> * ...
                     * ddc::type_seq_size_v<ValidIndexSetFirstDim>)>,
            ValidIndexSetFirstDim,
            ValidIndexSet...>;

public:
    using typename base_type::index_set;

    using base_type::rank;

private:
    using base_type::s_n_elements;
    std::array<ElementType, s_n_elements> m_data_alloc;
    using base_type::m_data;

public:
    KOKKOS_DEFAULTED_FUNCTION Tensor() = default;

    explicit KOKKOS_FUNCTION Tensor(ElementType fill_value)
    {
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data()[i] = fill_value;
        }
    }

    template <
            class... Params,
            class = std::enable_if_t<(std::is_convertible_v<Params, ElementType> && ...)>,
            class = std::enable_if_t<
                    sizeof...(Params) == base_type::size() && sizeof...(Params) != 1>>
    explicit KOKKOS_FUNCTION Tensor(Params... elements)
    {
        static_assert(
                rank() == 1,
                "Filling the tensor on initialisation is only permitted for 1D vector objects");
        m_data.m_data_alloc = std::array<ElementType, base_type::size()>({elements...});
    }

    template <class... Dims>
    explicit KOKKOS_FUNCTION Tensor(Coord<Dims...> coord) noexcept
    {
        static_assert(
                rank() == 1,
                "Filling the tensor on initialisation is only permitted for 1D vector objects");
        static_assert(
                std::is_same_v<VectorIndexSet<Dims...>, ddc::type_seq_element_t<0, index_set>>,
                "The coordinate must have the same memory layout to make a clean conversion.");
        m_data.m_data_alloc = coord.array();
    }

    template <class OTensorType, std::enable_if_t<is_tensor_type_v<OTensorType>, bool> = true>
    explicit KOKKOS_FUNCTION Tensor(const OTensorType& o_tensor) noexcept
    {
        static_assert(
                std::is_same_v<typename OTensorType::index_set, index_set>,
                "The coordinate must have the same memory layout to make a clean conversion.");
        for (std::size_t i(0); i < s_n_elements; ++i) {
            m_data()[i] = o_tensor.m_data()[i];
        }
    }

    KOKKOS_DEFAULTED_FUNCTION Tensor(Tensor const& o_tensor) = default;

    KOKKOS_DEFAULTED_FUNCTION Tensor(Tensor&& o_tensor) = default;

    KOKKOS_DEFAULTED_FUNCTION Tensor& operator=(Tensor const& other) = default;

    KOKKOS_DEFAULTED_FUNCTION Tensor& operator=(Tensor&& other) = default;
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

namespace detail {

template <class ElementType, class... ValidIndexSet>
inline constexpr bool enable_tensor_type<Tensor<ElementType, ValidIndexSet...>> = true;

}
```


