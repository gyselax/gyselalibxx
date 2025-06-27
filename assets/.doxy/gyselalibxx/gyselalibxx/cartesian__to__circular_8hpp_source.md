

# File cartesian\_to\_circular.hpp

[**File List**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**cartesian\_to\_circular.hpp**](cartesian__to__circular_8hpp.md)

[Go to the documentation of this file](cartesian__to__circular_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_aliases.hpp"
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
    using CoordJacobian = CoordArg;

    using X_cov = typename X::Dual;
    using Y_cov = typename Y::Dual;
    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;

private:
    Coord<X, Y> m_o_point;

public:
    explicit KOKKOS_FUNCTION CartesianToCircular(Coord<X, Y> o_point = Coord<X, Y>(0.0, 0.0))
        : m_o_point(o_point)
    {
    }

    KOKKOS_DEFAULTED_FUNCTION CartesianToCircular(CartesianToCircular const& other) = default;

    CartesianToCircular(CartesianToCircular&& x) = default;

    KOKKOS_DEFAULTED_FUNCTION ~CartesianToCircular() = default;

    CartesianToCircular& operator=(CartesianToCircular const& x) = default;

    CartesianToCircular& operator=(CartesianToCircular&& x) = default;

    KOKKOS_FUNCTION Coord<R, Theta> operator()(Coord<X, Y> const& coord) const
    {
        const double x = ddc::get<X>(coord) - ddc::get<X>(m_o_point);
        const double y = ddc::get<Y>(coord) - ddc::get<Y>(m_o_point);
        const double r = Kokkos::sqrt(x * x + y * y);
        const double theta_pi_to_pi(Kokkos::atan2(y, x));
        const double theta = theta_pi_to_pi * (theta_pi_to_pi >= 0)
                             + (theta_pi_to_pi + 2 * M_PI) * (theta_pi_to_pi < 0);
        return Coord<R, Theta>(r, theta);
    }

    KOKKOS_FUNCTION double jacobian(Coord<X, Y> const& coord)
    {
        const double x = ddc::get<X>(coord) - ddc::get<X>(m_o_point);
        const double y = ddc::get<Y>(coord) - ddc::get<Y>(m_o_point);
        return 1. / Kokkos::sqrt(x * x + y * y);
    }

    KOKKOS_FUNCTION DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X_cov, Y_cov>> jacobian_matrix(
            Coord<X, Y> const& coord) const

    {
        const double x = ddc::get<X>(coord) - ddc::get<X>(m_o_point);
        const double y = ddc::get<Y>(coord) - ddc::get<Y>(m_o_point);
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

        const double x = ddc::get<X>(coord) - ddc::get<X>(m_o_point);
        const double y = ddc::get<Y>(coord) - ddc::get<Y>(m_o_point);
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

    KOKKOS_INLINE_FUNCTION CircularToCartesian<R, Theta, X, Y> get_inverse_mapping() const
    {
        return CircularToCartesian<R, Theta, X, Y>(m_o_point);
    }
};

namespace mapping_detail {
template <class X, class Y, class R, class Theta, class ExecSpace>
struct MappingAccessibility<ExecSpace, CartesianToCircular<X, Y, R, Theta>> : std::true_type
{
};
} // namespace mapping_detail
```


