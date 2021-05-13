#pragma once

#include <array>
#include <ostream>
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


template <class ElementType>
struct accessor
{
    using offset_policy = accessor;
    using element_type = std::remove_cv_t<ElementType>;
    using reference = ElementType&;
    using pointer = ElementType*;

    MDSPAN_INLINE_FUNCTION
    constexpr accessor() noexcept = default;

    MDSPAN_TEMPLATE_REQUIRES(
            class OtherElementType,
            /* requires */ (_MDSPAN_TRAIT(std::is_convertible, OtherElementType, ElementType)))
    MDSPAN_INLINE_FUNCTION
    constexpr accessor(accessor<OtherElementType>&&) noexcept {}

    MDSPAN_TEMPLATE_REQUIRES(
            class OtherElementType,
            /* requires */ (_MDSPAN_TRAIT(std::is_convertible, OtherElementType, ElementType)))
    MDSPAN_INLINE_FUNCTION
    constexpr accessor(const accessor<OtherElementType>&) noexcept {}

    MDSPAN_INLINE_FUNCTION
    constexpr pointer offset(pointer p, ptrdiff_t i) const noexcept
    {
        return p + i;
    }

    MDSPAN_FORCE_INLINE_FUNCTION
    constexpr reference access(pointer p, ptrdiff_t i) const noexcept
    {
        return p[i];
    }

    MDSPAN_INLINE_FUNCTION
    constexpr pointer decay(pointer p) const noexcept
    {
        return p;
    }
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
    using type = std::experimental::basic_mdspan<
            ElementType,
            ExtentsND<N>,
            std::experimental::layout_right,
            detail::accessor<ElementType>>;
};

template <int N, class ElementType>
struct ViewNDMaker<N, ElementType, false>
{
    using type = std::experimental::
            basic_mdspan<ElementType, ExtentsND<N>, NCLayoutND<N>, detail::accessor<ElementType>>;
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


template <class>
struct IsContiguous;
template <ptrdiff_t S>
struct IsContiguous<std::experimental::layout_stride<S>>
{
    static constexpr bool val = (S == 1);
};
template <ptrdiff_t SH, ptrdiff_t SN, ptrdiff_t... ST>
struct IsContiguous<std::experimental::layout_stride<SH, SN, ST...>>
{
    static constexpr bool val = IsContiguous<std::experimental::layout_stride<SN, ST...>>::val;
};
template <>
struct IsContiguous<std::experimental::layout_right>
{
    static constexpr bool val = true;
};
template <class ElementType, class Extents, class LayoutPolicy, class AccessorPolicy>
struct IsContiguous<
        std::experimental::basic_mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>>
{
    static constexpr bool val = IsContiguous<LayoutPolicy>::val;
};

template <class C>
constexpr bool is_contiguous = IsContiguous<C>::val;

namespace details {

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
        std::experimental::basic_mdspan<ElementType, Extents, Layout, Accessor> const& s,
        std::index_sequence<I0, Is...>)
{
    if constexpr (sizeof...(Is) > 0) {
        os << '[';
        for (ptrdiff_t i0 = 0; i0 < s.extent(I0); ++i0) {
            stream_impl(
                    os,
                    std::experimental::subspan(s, i0, (Is, std::experimental::all)...),
                    std::make_index_sequence<sizeof...(Is)>());
        }
        os << ']';
    } else {
        os << '[';
        for (ptrdiff_t i0 = 0; i0 < s.extent(I0) - 1; ++i0) {
            os << s(i0) << ',';
        }
        os << s(s.extent(I0) - 1) << ']';
    }
    return os;
}

} // namespace details

/// Convenient function to dump a basic_mdspan, it recursively prints all dimensions.
/// Disclaimer: use with caution for large arrays
template <class ElementType, class Extents, class Layout, class Accessor>
std::ostream& operator<<(
        std::ostream& os,
        std::experimental::basic_mdspan<ElementType, Extents, Layout, Accessor> const& s)
{
    return details::stream_impl(os, s, std::make_index_sequence<Extents::rank()>());
}
