

# File cartesian\_to\_czarny.hpp

[**File List**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**cartesian\_to\_czarny.hpp**](cartesian__to__czarny_8hpp.md)

[Go to the documentation of this file](cartesian__to__czarny_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_aliases.hpp"
#include "view.hpp"

// Pre-declaration of analytical inverse
template <class R, class Theta, class X, class Y>
class CzarnyToCartesian;

template <class X, class Y, class R, class Theta>
class CartesianToCzarny
{
public:
    using cartesian_tag_x = X;
    using cartesian_tag_y = Y;
    using curvilinear_tag_r = R;
    using curvilinear_tag_theta = Theta;

    using CoordArg = Coord<X, Y>;
    using CoordResult = Coord<R, Theta>;

private:
    double m_epsilon;
    double m_e;
    double m_x0;
    double m_y0;

public:
    explicit KOKKOS_FUNCTION CartesianToCzarny(
            double epsilon,
            double e,
            double x0 = 0.0,
            double y0 = 0.0)
        : m_epsilon(epsilon)
        , m_e(e)
        , m_x0(x0)
        , m_y0(y0)
    {
    }

    KOKKOS_DEFAULTED_FUNCTION CartesianToCzarny(CartesianToCzarny const& other) = default;

    CartesianToCzarny(CartesianToCzarny&& x) = default;

    KOKKOS_DEFAULTED_FUNCTION ~CartesianToCzarny() = default;

    CartesianToCzarny& operator=(CartesianToCzarny const& x) = default;

    CartesianToCzarny& operator=(CartesianToCzarny&& x) = default;

    KOKKOS_FUNCTION double epsilon() const
    {
        return m_epsilon;
    }

    KOKKOS_FUNCTION double e() const
    {
        return m_e;
    }

    KOKKOS_FUNCTION Coord<R, Theta> operator()(Coord<X, Y> const& coord) const
    {
        const double x = ddc::get<X>(coord) - m_x0;
        const double y = ddc::get<Y>(coord) - m_y0;
        const double ex = 1. + m_epsilon * x;
        const double ex2 = (m_epsilon * x * x - 2. * x - m_epsilon);
        const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
        const double xi = Kokkos::sqrt(xi2);
        const double r = Kokkos::sqrt(y * y * ex * ex / (m_e * m_e * xi2) + ex2 * ex2 * 0.25);
        double theta = Kokkos::atan2(2. * y * ex, (m_e * xi * ex2));
        if (theta < 0) {
            theta = 2 * M_PI + theta;
        }
        return Coord<R, Theta>(r, theta);
    }

    KOKKOS_INLINE_FUNCTION CzarnyToCartesian<R, Theta, X, Y> get_inverse_mapping() const
    {
        return CzarnyToCartesian<R, Theta, X, Y>(m_epsilon, m_e, m_x0, m_y0);
    }
};

namespace mapping_detail {
template <class X, class Y, class R, class Theta, class ExecSpace>
struct MappingAccessibility<ExecSpace, CartesianToCzarny<X, Y, R, Theta>> : std::true_type
{
};
} // namespace mapping_detail
```


