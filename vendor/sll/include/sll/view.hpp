#pragma once

#include <array>
#include <ostream>
#include <utility>

#include <experimental/mdspan>

namespace detail {

template <std::size_t N, class ElementType, bool CONTIGUOUS = true>
struct ViewNDMaker;

template <std::size_t N, class ElementType>
struct ViewNDMaker<N, ElementType, true>
{
    using type = std::experimental::
            mdspan<ElementType, std::experimental::dextents<N>, std::experimental::layout_right>;
};

template <std::size_t N, class ElementType>
struct ViewNDMaker<N, ElementType, false>
{
    using type = std::experimental::
            mdspan<ElementType, std::experimental::dextents<N>, std::experimental::layout_stride>;
};

/// Note: We use the comma operator to fill the input parameters
///
/// If Is=[1, 2], `subspan(s, i0, (Is, all)...)` will be expanded as
/// `subspan(s, i0, (1, all), (2, all))` which is equivalent to
/// `subspan(s, i0, all, all)`
template <
        class ElementType,
        class Extents,
        class Layout,
        class Accessor,
        std::size_t I0,
        std::size_t... Is>
std::ostream& stream_impl(
        std::ostream& os,
        std::experimental::mdspan<ElementType, Extents, Layout, Accessor> const& s,
        std::index_sequence<I0, Is...>)
{
    if constexpr (sizeof...(Is) > 0) {
        os << '[';
        for (std::size_t i0 = 0; i0 < s.extent(I0); ++i0) {
            stream_impl(
                    os,
                    std::experimental::
                            submdspan(s, i0, ((void)Is, std::experimental::full_extent)...),
                    std::make_index_sequence<sizeof...(Is)>());
        }
        os << ']';
    } else {
        os << '[';
        for (std::size_t i0 = 0; i0 < s.extent(I0) - 1; ++i0) {
            os << s(i0) << ',';
        }
        os << s(s.extent(I0) - 1) << ']';
    }
    return os;
}


// template <class>
// struct IsContiguous;
// template <>
// struct IsContiguous<std::experimental::layout_stride>
// {
//     static constexpr bool val = false;
// };
// template <>
// struct IsContiguous<std::experimental::layout_right>
// {
//     static constexpr bool val = true;
// };
// template <class ElementType, class Extents, class LayoutPolicy, class AccessorPolicy>
// struct IsContiguous<std::experimental::mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>>
// {
//     static constexpr bool val = IsContiguous<LayoutPolicy>::val;
// };

} // namespace detail


template <std::size_t N, class ElementType>
using SpanND = std::experimental::mdspan<ElementType, std::experimental::dextents<N>>;

template <std::size_t N, class ElementType>
using ViewND = SpanND<N, ElementType const>;

template <class ElementType>
using Span1D = SpanND<1, ElementType>;

template <class ElementType>
using Span2D = SpanND<2, ElementType>;

template <class ElementType>
using View1D = ViewND<1, ElementType>;

template <class ElementType>
using View2D = ViewND<2, ElementType>;

using DSpan1D = Span1D<double>;

using DSpan2D = Span2D<double>;

using DView1D = View1D<double>;

using DView2D = View2D<double>;

// template <class C>
// constexpr bool is_contiguous_v = detail::IsContiguous<C>::val;

/// Convenient function to dump a mdspan, it recursively prints all dimensions.
/// Disclaimer: use with caution for large arrays
template <class ElementType, class Extents, class Layout, class Accessor>
std::ostream& operator<<(
        std::ostream& os,
        std::experimental::mdspan<ElementType, Extents, Layout, Accessor> const& s)
{
    return detail::stream_impl(os, s, std::make_index_sequence<Extents::rank()>());
}
