

# File cartesian\_to\_circular.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**cartesian\_to\_circular.hpp**](cartesian__to__circular_8hpp.md)

[Go to the documentation of this file](cartesian__to__circular_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "tensor.hpp"
#include "view.hpp"

// Pre-declaration of analytical inverse
template <class R, class Theta, class X, class Y>
class CircularToCartesian;

template <class X, class Y, class R, class Theta>
class CartesianToCircular
{
public:
    using cartesian_tag_x = X;
    using cartesian_tag_y = Y;
    using curvilinear_tag_r = R;
    using curvilinear_tag_theta = Theta;

    using CoordArg = Coord<X, Y>;
    using CoordResult = Coord<R, Theta>;

    using X_cov = typename X::Dual;
    using Y_cov = typename Y::Dual;
    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;

public:
    CartesianToCircular() = default;

    KOKKOS_FUNCTION CartesianToCircular(CartesianToCircular const& other) {}

    CartesianToCircular(CartesianToCircular&& x) = default;

    ~CartesianToCircular() = default;

    CartesianToCircular& operator=(CartesianToCircular const& x) = default;

    CartesianToCircular& operator=(CartesianToCircular&& x) = default;

    KOKKOS_FUNCTION Coord<R, Theta> operator()(Coord<X, Y> const& coord) const
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        const double r = Kokkos::sqrt(x * x + y * y);
        const double theta = Kokkos::atan2(y, x);
        return Coord<R, Theta>(r, theta);
    }

    KOKKOS_FUNCTION double jacobian(Coord<X, Y> const& coord)
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        return 2. * (x * x - y * y) / Kokkos::pow(x * x + y * y, 1.5);
    }

    KOKKOS_FUNCTION DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix(
            Coord<X, Y> const& coord) const

    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X_cov, Y_cov>> jacobian_matrix;
        ddcHelper::get<R, X_cov>(jacobian_matrix) = 2 * x / Kokkos::pow(x * x + y * y, 0.5);
        ddcHelper::get<R, Y_cov>(jacobian_matrix) = 2 * y / Kokkos::pow(x * x + y * y, 0.5);
        ddcHelper::get<Theta, X_cov>(jacobian_matrix) = -y / Kokkos::pow(x * x + y * y, 2.);
        ddcHelper::get<Theta, Y_cov>(jacobian_matrix) = x / Kokkos::pow(x * x + y * y, 2.);
        return jacobian_matrix;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double jacobian_component(Coord<X, Y> const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<X, Y>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<R_cov, Theta_cov>>);

        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (1,1), i.e dx/dr
            return 2 * x / Kokkos::pow(x * x + y * y, 0.5);
        } else if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, Theta_cov>) {
            // Component (1,2), i.e dx/dtheta
            return 2 * y / Kokkos::pow(x * x + y * y, 0.5);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (2,1), i.e dy/dr
            return -y / Kokkos::pow(x * x + y * y, 2.);
        } else {
            // Component (2,2), i.e dy/dtheta
            return x / Kokkos::pow(x * x + y * y, 2.);
        }
    }

    CircularToCartesian<R, Theta, X, Y> get_inverse_mapping() const
    {
        return CircularToCartesian<R, Theta, X, Y>();
    }
};

namespace mapping_detail {
template <class X, class Y, class R, class Theta, class ExecSpace>
struct MappingAccessibility<ExecSpace, CartesianToCircular<X, Y, R, Theta>> : std::true_type
{
};
} // namespace mapping_detail
```


