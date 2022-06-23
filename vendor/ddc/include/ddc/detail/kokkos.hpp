#pragma once

#include <experimental/mdspan>

#include <Kokkos_Core.hpp>

#include "macros.hpp"

namespace detail {

template <class KokkosLP>
struct mdspan_layout;

template <class KokkosLP>
using mdspan_layout_t = typename mdspan_layout<KokkosLP>::type;

template <>
struct mdspan_layout<Kokkos::LayoutLeft>
{
    using type = std::experimental::layout_left;
};

template <>
struct mdspan_layout<Kokkos::LayoutRight>
{
    using type = std::experimental::layout_right;
};

template <>
struct mdspan_layout<Kokkos::LayoutStride>
{
    using type = std::experimental::layout_stride;
};


template <class mdspanLP>
struct kokkos_layout;

template <class mdspanLP>
using kokkos_layout_t = typename kokkos_layout<mdspanLP>::type;

template <>
struct kokkos_layout<std::experimental::layout_left>
{
    using type = Kokkos::LayoutLeft;
};

template <>
struct kokkos_layout<std::experimental::layout_right>
{
    using type = Kokkos::LayoutRight;
};

template <>
struct kokkos_layout<std::experimental::layout_stride>
{
    using type = Kokkos::LayoutStride;
};


template <std::size_t... Is>
Kokkos::LayoutStride make_layout_stride(
        std::array<std::size_t, sizeof...(Is)> const& interleaved_extents_strides,
        std::index_sequence<Is...>)
{
    return Kokkos::LayoutStride(interleaved_extents_strides[Is]...);
}

template <class EP, class MP, std::size_t... Is>
kokkos_layout_t<typename MP::layout_type> build_kokkos_layout(
        EP const& ep,
        MP const& mapping,
        std::index_sequence<Is...>)
{
    DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
    using kokkos_layout_type = kokkos_layout_t<typename MP::layout_type>;
    if constexpr (std::is_same_v<kokkos_layout_type, Kokkos::LayoutStride>) {
        std::array<std::size_t, sizeof...(Is) * 2> storage;
        std::experimental::mdspan<
                std::size_t,
                std::experimental::extents<sizeof...(Is), sizeof...(Is)>,
                std::experimental::layout_right>
                interleaved_extents_strides(storage.data());
        ((interleaved_extents_strides(Is, 0) = ep.extent(Is),
          interleaved_extents_strides(Is, 1) = mapping.stride(Is)),
         ...);
        return make_layout_stride(storage, std::make_index_sequence<sizeof...(Is) * 2> {});
    } else {
        return kokkos_layout_type(ep.extent(Is)...);
    }
    DDC_IF_NVCC_THEN_POP
}

/// Recursively add a pointer
template <class ET, std::size_t N>
struct mdspan_to_kokkos_element_type : mdspan_to_kokkos_element_type<std::add_pointer_t<ET>, N - 1>
{
};

template <class ET>
struct mdspan_to_kokkos_element_type<ET, 0>
{
    using type = ET;
};

template <class ET, std::size_t N>
using mdspan_to_kokkos_element_type_t = typename mdspan_to_kokkos_element_type<ET, N>::type;

template <class T, std::size_t N>
struct final_type
{
    using type = T;
    static constexpr std::size_t rank = N;
};

/// Recursively remove a pointer
template <class ET, std::size_t N>
struct kokkos_to_mdspan_element_type
    : std::conditional_t<
              std::is_pointer_v<std::decay_t<ET>>,
              kokkos_to_mdspan_element_type<std::remove_pointer_t<std::decay_t<ET>>, N + 1>,
              final_type<ET, N>>
{
};

template <class ET>
using kokkos_to_mdspan_element_type_t = typename kokkos_to_mdspan_element_type<ET, 0>::type;

template <class ET>
constexpr inline std::size_t kokkos_to_mdspan_element_type_rank
        = kokkos_to_mdspan_element_type<ET, 0>::rank;

template <class DataType, class... Properties, std::size_t... Is>
auto build_mdspan(Kokkos::View<DataType, Properties...> const view, std::index_sequence<Is...>)
{
    DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
    using element_type = kokkos_to_mdspan_element_type_t<DataType>;
    using extents_type = std::experimental::dextents<Kokkos::View<DataType, Properties...>::rank>;
    using layout_type
            = mdspan_layout_t<typename Kokkos::View<DataType, Properties...>::array_layout>;
    using mapping_type = typename layout_type::template mapping<extents_type>;
    extents_type exts(view.extent(Is)...);
    if constexpr (std::is_same_v<layout_type, std::experimental::layout_stride>) {
        return std::experimental::mdspan(view.data(), mapping_type(exts, {view.stride(Is)...}));
    } else {
        return std::experimental::
                mdspan<element_type, extents_type, layout_type>(view.data(), mapping_type(exts));
    }
    DDC_IF_NVCC_THEN_POP
}

} // namespace detail
