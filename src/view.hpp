#pragma once

#include <experimental/mdspan>

namespace detail {

template <int N, typename, typename = std::enable_if<(N >= 0)>>
struct ExtentAddND;

template <typename Extents>
struct ExtentAddND<0, Extents> {
    using type = Extents;
};

template <int N, ptrdiff_t... Extents>
struct ExtentAddND<N, std::experimental::extents<Extents...>, std::enable_if<(N > 0)>> {
    using type = typename ExtentAddND<N - 1, std::experimental::extents<std::experimental::dynamic_extent, Extents...>>::type;
};

template <int N>
using ExtentND = typename ExtentAddND<N, std::experimental::extents<>>::type;

template <int N, typename LayoutPolicy, typename = std::enable_if<(N >= 0)>>
struct LayoutAddND;

template <typename LayoutPolicy>
struct LayoutAddND<0, LayoutPolicy> {
    using type = LayoutPolicy;
};

template <int N, ptrdiff_t... Strides>
struct LayoutAddND<N, std::experimental::layout_stride<Strides...>, std::enable_if<(N > 0)>> {
    using type = typename LayoutAddND<N - 1, std::experimental::layout_stride<std::experimental::dynamic_extent, Strides...>>::type;
};

template <int N>
using LayoutND = typename LayoutAddND<N - 1, std::experimental::layout_stride<1>>::type;

}

template <int N>
using ViewND = std::experimental::basic_mdspan<double, detail::ExtentND<N>, detail::LayoutND<N>>;
