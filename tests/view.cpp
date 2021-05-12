#include <span>

#include <gtest/gtest.h>

#include "view.h"

using namespace std;
using namespace std::experimental;


template <class ElementType>
struct accessor
{
    using offset_policy = accessor;
    using element_type = std::remove_cv_t<ElementType>;
    using reference = ElementType &;
    using pointer = ElementType *;

    MDSPAN_INLINE_FUNCTION
    constexpr explicit accessor() noexcept = default;

    MDSPAN_INLINE_FUNCTION
    constexpr accessor(accessor<element_type> const&) noexcept {}

    MDSPAN_INLINE_FUNCTION
    constexpr accessor(accessor<element_type>&&) noexcept {}

    MDSPAN_INLINE_FUNCTION
    constexpr accessor(accessor<element_type const> const&) noexcept {}

    MDSPAN_INLINE_FUNCTION
    constexpr accessor(accessor<element_type const>&&) noexcept {}

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

TEST(View1DTest, Constructor)
{
    std::array<double, 10> x = {0};
    basic_mdspan<double, extents<std::experimental::dynamic_extent>, layout_right, accessor<double>> xv(x.data(), x.size());
    basic_mdspan<double, extents<std::experimental::dynamic_extent>, layout_right, accessor<double>> xv_(xv);
    basic_mdspan<double, extents<std::experimental::dynamic_extent>, layout_right, accessor<double const>> xcv(xv);
    basic_mdspan<double, extents<std::experimental::dynamic_extent>, layout_right, accessor<double const>> xcv_(xcv);
    std::array<double, 10> const cx = {0};
    basic_mdspan<double, extents<std::experimental::dynamic_extent>, layout_right, accessor<double const>> cxcv(cx.data(), cx.size());
    basic_mdspan<double, extents<std::experimental::dynamic_extent>, layout_right, accessor<double const>> cxcv_(cxcv);
    std::array<double const, 10> const cx2 = {0};
    basic_mdspan<double, extents<std::experimental::dynamic_extent>, layout_right, accessor<double const>> cxcv2(cx2.data(), cx2.size());
    basic_mdspan<double, extents<std::experimental::dynamic_extent>, layout_right, accessor<double const>> cxcv2_(cxcv2);
}
