

# File cartesian\_to\_cylindrical.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**cartesian\_to\_cylindrical.hpp**](cartesian__to__cylindrical_8hpp.md)

[Go to the documentation of this file](cartesian__to__cylindrical_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "tensor.hpp"
#include "view.hpp"

// Pre-declaration of analytical inverse
template <class R, class Z, class Zeta, class X, class Y>
class CylindricalToCartesian;

template <class X, class Y, class Z, class R, class Zeta>
class CartesianToCylindrical
{
public:
    using cartesian_tag_x = X;
    using cartesian_tag_y = Y;
    using cartesian_tag_z = Z;
    using cylindrical_tag_R = R;
    using cylindrical_tag_Z = Z;
    using cylindrical_tag_Zeta = Zeta;

    using CoordArg = Coord<X, Y, Z>;
    using CoordResult = Coord<R, Z, Zeta>;

    using X_cov = typename X::Dual;
    using Y_cov = typename Y::Dual;
    using Z_cov = typename Z::Dual;
    using R_cov = typename R::Dual;
    using Zeta_cov = typename Zeta::Dual;

public:
    CartesianToCylindrical() = default;

    KOKKOS_FUNCTION CartesianToCylindrical(CartesianToCylindrical const& other) {}

    CartesianToCylindrical(CartesianToCylindrical&& x) = default;

    ~CartesianToCylindrical() = default;

    CartesianToCylindrical& operator=(CartesianToCylindrical const& x) = default;

    CartesianToCylindrical& operator=(CartesianToCylindrical&& x) = default;

    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        Coord<R> r(Kokkos::sqrt(x * x + y * y));
        const double zeta_pi_to_pi(Kokkos::atan2(y, x));
        Coord<Zeta> zeta(zeta_pi_to_pi + 2 * M_PI * (zeta_pi_to_pi < 0));

        return CoordResult(r, ddc::select<Z>(coord), zeta);
    }

    KOKKOS_FUNCTION double jacobian(CoordArg const& coord)
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        return -1. / Kokkos::sqrt(x * x + y * y);
    }

    KOKKOS_FUNCTION DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<X_cov, Y_cov, Z_cov>>
    jacobian_matrix(CoordArg const& coord) const
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);

        DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<X_cov, Y_cov, Z_cov>> matrix;
        ddcHelper::get<R, X_cov>(matrix) = x / Kokkos::sqrt(x * x + y * y);
        ddcHelper::get<R, Y_cov>(matrix) = y / Kokkos::sqrt(x * x + y * y);
        ddcHelper::get<R, Z_cov>(matrix) = 0;
        ddcHelper::get<Z, X_cov>(matrix) = 0;
        ddcHelper::get<Z, Y_cov>(matrix) = 0;
        ddcHelper::get<Z, Z_cov>(matrix) = 1;
        ddcHelper::get<Zeta, X_cov>(matrix) = -y / (x * x + y * y);
        ddcHelper::get<Zeta, Y_cov>(matrix) = x / (x * x + y * y);
        ddcHelper::get<Zeta, Z_cov>(matrix) = 0;

        return matrix;
    }


    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<X, Y>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<R_cov, Zeta_cov>>);
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<R, Z, Zeta>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<X_cov, Y_cov, Z_cov>>);

        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (1,1), i.e dx/dr
            return x / Kokkos::sqrt(x * x + y * y);
        } else if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, Zeta_cov>) {
            // Component (1,2), i.e dx/dzeta
            return y / Kokkos::sqrt(x * x + y * y);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (2,1), i.e dy/dr
            return -y / (x * x + y * y);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, Zeta_cov>) {
            // Component (2,2), i.e dy/dzeta
            return x / (x * x + y * y);
        } else if constexpr (std::is_same_v<IndexTag1, Z> && std::is_same_v<IndexTag2, Z_cov>) {
            return 1;
        } else {
            return 0;
        }
    }

    CylindricalToCartesian<R, Z, Zeta, X, Y> get_inverse_mapping() const
    {
        return CylindricalToCartesian<R, Z, Zeta, X, Y>();
    }
};
```


