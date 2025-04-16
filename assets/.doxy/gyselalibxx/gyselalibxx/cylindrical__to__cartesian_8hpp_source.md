

# File cylindrical\_to\_cartesian.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**cylindrical\_to\_cartesian.hpp**](cylindrical__to__cartesian_8hpp.md)

[Go to the documentation of this file](cylindrical__to__cartesian_8hpp.md)


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
template <class X, class Y, class Z, class R, class Zeta>
class CartesianToCylindrical;

template <class R, class Z, class Zeta, class X, class Y>
class CylindricalToCartesian
{
public:
    using cartesian_tag_x = X;
    using cartesian_tag_y = Y;
    using cartesian_tag_z = Z;
    using cylindrical_tag_R = R;
    using cylindrical_tag_Z = Z;
    using cylindrical_tag_Zeta = Zeta;

    using CoordArg = Coord<R, Z, Zeta>;
    using CoordResult = Coord<X, Y, Z>;

    using X_cov = typename X::Dual;
    using Y_cov = typename Y::Dual;
    using Z_cov = typename Z::Dual;
    using R_cov = typename R::Dual;
    using Zeta_cov = typename Zeta::Dual;

public:
    CylindricalToCartesian() = default;

    KOKKOS_FUNCTION CylindricalToCartesian(CylindricalToCartesian const& other) {}

    CylindricalToCartesian(CylindricalToCartesian&& x) = default;

    ~CylindricalToCartesian() = default;

    CylindricalToCartesian& operator=(CylindricalToCartesian const& x) = default;

    CylindricalToCartesian& operator=(CylindricalToCartesian&& x) = default;

    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double zeta = ddc::get<Zeta>(coord);
        Coord<X> x(r * Kokkos::cos(zeta));
        Coord<Y> y(r * Kokkos::sin(zeta));
        return CoordResult(x, y, ddc::select<Z>(coord));
    }

    KOKKOS_FUNCTION double jacobian(CoordArg const& coord) const
    {
        double r = ddc::get<R>(coord);
        return -r;
    }

    KOKKOS_FUNCTION DTensor<VectorIndexSet<X, Y, Z>, VectorIndexSet<R_cov, Z_cov, Zeta_cov>>
    jacobian_matrix(CoordArg const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double zeta = ddc::get<Zeta>(coord);
        DTensor<VectorIndexSet<X, Y, Z>, VectorIndexSet<R_cov, Z_cov, Zeta_cov>> jacobian_matrix;
        ddcHelper::get<X, R_cov>(jacobian_matrix) = Kokkos::cos(zeta);
        ddcHelper::get<X, Z_cov>(jacobian_matrix) = 0;
        ddcHelper::get<X, Zeta_cov>(jacobian_matrix) = -r * Kokkos::sin(zeta);
        ddcHelper::get<Y, R_cov>(jacobian_matrix) = Kokkos::sin(zeta);
        ddcHelper::get<Y, Z_cov>(jacobian_matrix) = 0;
        ddcHelper::get<Y, Zeta_cov>(jacobian_matrix) = r * Kokkos::cos(zeta);
        ddcHelper::get<Z, R_cov>(jacobian_matrix) = 0;
        ddcHelper::get<Z, Z_cov>(jacobian_matrix) = 1;
        ddcHelper::get<Z, Zeta_cov>(jacobian_matrix) = 0;
        return jacobian_matrix;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<X, Y, Z>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<R_cov, Z_cov, Zeta_cov>>);

        const double zeta = ddc::get<Zeta>(coord);

        if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, R_cov>) {
            //Compute the (1,1) coefficient of the Jacobian matrix, i.e J^x_r.
            return Kokkos::cos(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, Zeta_cov>) {
            //Compute the (1,2) coefficient of the Jacobian matrix, i.e J^x_theta.
            const double r = ddc::get<R>(coord);
            return -r * Kokkos::sin(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, R_cov>) {
            //Compute the (2,1) coefficient of the Jacobian matrix, i.e J^y_r.
            return Kokkos::sin(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, Zeta_cov>) {
            //Compute the (2,2) coefficient of the Jacobian matrix, i.e J^y_theta.
            const double r = ddc::get<R>(coord);
            return r * Kokkos::cos(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, Z> && std::is_same_v<IndexTag2, Z_cov>) {
            return 1;
        } else {
            return 0;
        }
    }


    KOKKOS_FUNCTION DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<X_cov, Y_cov, Z_cov>>
    inv_jacobian_matrix(CoordArg const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double zeta = ddc::get<Zeta>(coord);
        assert(fabs(r) >= 1e-15);

        DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<X_cov, Y_cov, Z_cov>> matrix(0);
        ddcHelper::get<R, X_cov>(matrix) = Kokkos::cos(zeta);
        ddcHelper::get<R, Y_cov>(matrix) = Kokkos::sin(zeta);
        ddcHelper::get<Zeta, X_cov>(matrix) = -1 / r * Kokkos::sin(zeta);
        ddcHelper::get<Zeta, Y_cov>(matrix) = 1 / r * Kokkos::cos(zeta);
        ddcHelper::get<Z, Z_cov>(matrix) = 1;
        return matrix;
    }


    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double inv_jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<R, Z, Zeta>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<X_cov, Y_cov, Z_cov>>);

        const double zeta = ddc::get<Zeta>(coord);

        if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, X_cov>) {
            //Compute the (1,1) coefficient of the inverse Jacobian matrix.
            return Kokkos::cos(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, Y_cov>) {
            //Compute the (1,2) coefficient of the inverse Jacobian matrix.
            return Kokkos::sin(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, Z> && std::is_same_v<IndexTag2, Z_cov>) {
            return 1;
        } else if constexpr (std::is_same_v<IndexTag1, Z> || std::is_same_v<IndexTag2, Z_cov>) {
            return 0;
        } else {
            const double r = ddc::get<R>(coord);
            assert(fabs(r) >= 1e-15);
            if constexpr (std::is_same_v<IndexTag2, X_cov>) {
                //Compute the (2,1) coefficient of the inverse Jacobian matrix.
                return -1 / r * Kokkos::sin(zeta);
            } else {
                //Compute the (2,2) coefficient of the inverse Jacobian matrix.
                return 1 / r * Kokkos::cos(zeta);
            }
        }
    }


    CartesianToCylindrical<X, Y, Z, R, Zeta> get_inverse_mapping() const
    {
        return CartesianToCylindrical<X, Y, Z, R, Zeta>();
    }
};

namespace mapping_detail {
template <class X, class Y, class Z, class R, class Zeta, class ExecSpace>
struct MappingAccessibility<ExecSpace, CylindricalToCartesian<R, Z, Zeta, X, Y>> : std::true_type
{
};

} // namespace mapping_detail
```


