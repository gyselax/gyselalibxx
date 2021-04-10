#pragma once

#include <utility>

#include <experimental/mdspan>

namespace detail {

template <size_t>
constexpr auto dynamic_extent = std::experimental::dynamic_extent;

template <class S>
struct ExtentNDBuilder;

template <std::size_t... S>
struct ExtentNDBuilder<std::index_sequence<S...>>
{
    using type = std::experimental::extents<dynamic_extent<S>...>;
};

} // namespace detail


template <size_t N>
using ExtentsND = typename detail::ExtentNDBuilder<std::make_index_sequence<N>>::type;

using Extents1D = ExtentsND<1>;

using Extents2D = ExtentsND<2>;


namespace detail {

template <int N, typename LayoutPolicy, typename = std::enable_if<(N >= 0)>>
struct LayoutAddND;

template <typename LayoutPolicy>
struct LayoutAddND<0, LayoutPolicy>
{
    using type = LayoutPolicy;
};

template <int N, ptrdiff_t... Strides>
struct LayoutAddND<N, std::experimental::layout_stride<Strides...>, std::enable_if<(N > 0)>>
{
    using type = typename LayoutAddND<
            N - 1,
            std::experimental::layout_stride<std::experimental::dynamic_extent, Strides...>>::type;
};

} // namespace detail


template <int N>
using LayoutND = typename detail::LayoutAddND<N - 1, std::experimental::layout_stride<1>>::type;

template <int N>
using NCLayoutND = typename detail::LayoutAddND<N, std::experimental::layout_stride<>>::type;


namespace detail {

template <int N, class ElementType, bool CONTIGUOUS = true>
struct ViewNDMaker;

template <int N, class ElementType>
struct ViewNDMaker<N, ElementType, true>
{
    using type = std::experimental::basic_mdspan<ElementType, ExtentsND<N>>;
};

template <int N, class ElementType>
struct ViewNDMaker<N, ElementType, false>
{
    using type = std::experimental::basic_mdspan<ElementType, ExtentsND<N>, NCLayoutND<N>>;
};

} // namespace detail


template <int N, class ElementType, bool CONTIGUOUS = true>
using ViewND = typename detail::ViewNDMaker<N, ElementType, CONTIGUOUS>::type;

template <class ElementType, bool CONTIGUOUS = true>
using View1D = ViewND<1, ElementType, CONTIGUOUS>;

using mdspan_1d = View1D<double>;

template <class ElementType, bool CONTIGUOUS = true>
using View2D = ViewND<2, ElementType, CONTIGUOUS>;

using mdspan_2d = View2D<double>;
