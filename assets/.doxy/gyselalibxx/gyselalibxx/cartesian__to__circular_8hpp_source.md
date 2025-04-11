

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
        const double theta_pi_to_pi(Kokkos::atan2(y, x));
        const double theta = theta_pi_to_pi * (theta_pi_to_pi >= 0)
                             + (theta_pi_to_pi + 2 * M_PI) * (theta_pi_to_pi < 0);
        return Coord<R, Theta>(r, theta);
    }

    KOKKOS_FUNCTION double jacobian(Coord<X, Y> const& coord)
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        return 1. / Kokkos::sqrt(x * x + y * y);
    }

    KOKKOS_FUNCTION DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X_cov, Y_cov>> jacobian_matrix(
            Coord<X, Y> const& coord) const

    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X_cov, Y_cov>> jacobian_matrix;
        const double r2 = x * x + y * y;
        const double r = Kokkos::sqrt(r2);
        ddcHelper::get<R, X_cov>(jacobian_matrix) = x / r;
        ddcHelper::get<R, Y_cov>(jacobian_matrix) = y / r;
        ddcHelper::get<Theta, X_cov>(jacobian_matrix) = -y / r2;
        ddcHelper::get<Theta, Y_cov>(jacobian_matrix) = x / r2;
        return jacobian_matrix;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double jacobian_component(Coord<X, Y> const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<R, Theta>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<X_cov, Y_cov>>);

        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, X_cov>) {
            // Component (1,1), i.e dr/dx
            return x / Kokkos::pow(x * x + y * y, 0.5);
        } else if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, Y_cov>) {
            // Component (1,2), i.e dr/dy
            return y / Kokkos::pow(x * x + y * y, 0.5);
        } else if constexpr (std::is_same_v<IndexTag1, Theta> && std::is_same_v<IndexTag2, X_cov>) {
            // Component (2,1), i.e dtheta/dy
            return -y / (x * x + y * y);
        } else {
            // Component (2,2), i.e dtheta/dy
            return x / (x * x + y * y);
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


