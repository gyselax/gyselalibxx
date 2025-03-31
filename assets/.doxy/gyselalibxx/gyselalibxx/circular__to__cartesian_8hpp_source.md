

# File circular\_to\_cartesian.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**circular\_to\_cartesian.hpp**](circular__to__cartesian_8hpp.md)

[Go to the documentation of this file](circular__to__cartesian_8hpp.md)


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
template <class X, class Y, class R, class Theta>
class CartesianToCircular;

template <class R, class Theta, class X, class Y>
class CircularToCartesian
{
public:
    using cartesian_tag_x = X;
    using cartesian_tag_y = Y;
    using curvilinear_tag_r = R;
    using curvilinear_tag_theta = Theta;

    using CoordArg = Coord<R, Theta>;
    using CoordResult = Coord<X, Y>;

    using X_cov = typename X::Dual;
    using Y_cov = typename Y::Dual;
    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;

public:
    CircularToCartesian() = default;

    KOKKOS_FUNCTION CircularToCartesian(CircularToCartesian const& other) {}

    CircularToCartesian(CircularToCartesian&& x) = default;

    ~CircularToCartesian() = default;

    CircularToCartesian& operator=(CircularToCartesian const& x) = default;

    CircularToCartesian& operator=(CircularToCartesian&& x) = default;

    KOKKOS_FUNCTION Coord<X, Y> operator()(Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        const double x = r * Kokkos::cos(theta);
        const double y = r * Kokkos::sin(theta);
        return Coord<X, Y>(x, y);
    }

    KOKKOS_FUNCTION double jacobian(Coord<R, Theta> const& coord) const
    {
        double r = ddc::get<R>(coord);
        return r;
    }

    KOKKOS_FUNCTION DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix(
            Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix;
        ddcHelper::get<X, R_cov>(jacobian_matrix) = Kokkos::cos(theta);
        ddcHelper::get<X, Theta_cov>(jacobian_matrix) = -r * Kokkos::sin(theta);
        ddcHelper::get<Y, R_cov>(jacobian_matrix) = Kokkos::sin(theta);
        ddcHelper::get<Y, Theta_cov>(jacobian_matrix) = r * Kokkos::cos(theta);
        return jacobian_matrix;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double jacobian_component(Coord<R, Theta> const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<X, Y>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<R_cov, Theta_cov>>);

        const double theta = ddc::get<Theta>(coord);

        if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, R_cov>) {
            //Compute the (1,1) coefficient of the Jacobian matrix, i.e J^x_r.
            return Kokkos::cos(theta);
        } else if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, Theta_cov>) {
            //Compute the (1,2) coefficient of the Jacobian matrix, i.e J^x_theta.
            const double r = ddc::get<R>(coord);
            return -r * Kokkos::sin(theta);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, R_cov>) {
            //Compute the (2,1) coefficient of the Jacobian matrix, i.e J^y_r.
            return Kokkos::sin(theta);
        } else {
            //Compute the (2,2) coefficient of the Jacobian matrix, i.e J^y_theta.
            const double r = ddc::get<R>(coord);
            return r * Kokkos::cos(theta);
        }
    }


    KOKKOS_FUNCTION DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X_cov, Y_cov>>
    inv_jacobian_matrix(Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        assert(fabs(r) >= 1e-15);

        DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X_cov, Y_cov>> matrix;
        ddcHelper::get<R, X_cov>(matrix) = Kokkos::cos(theta);
        ddcHelper::get<R, Y_cov>(matrix) = Kokkos::sin(theta);
        ddcHelper::get<Theta, X_cov>(matrix) = -1 / r * Kokkos::sin(theta);
        ddcHelper::get<Theta, Y_cov>(matrix) = 1 / r * Kokkos::cos(theta);
        return matrix;
    }


    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double inv_jacobian_component(Coord<R, Theta> const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<R, Theta>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<X_cov, Y_cov>>);

        const double theta = ddc::get<Theta>(coord);

        if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, X_cov>) {
            //Compute the (1,1) coefficient of the inverse Jacobian matrix.
            return Kokkos::cos(theta);
        } else if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, Y_cov>) {
            //Compute the (1,2) coefficient of the inverse Jacobian matrix.
            return Kokkos::sin(theta);
        } else {
            const double r = ddc::get<R>(coord);
            assert(fabs(r) >= 1e-15);
            if constexpr (std::is_same_v<IndexTag2, X_cov>) {
                //Compute the (2,1) coefficient of the inverse Jacobian matrix.
                return -1 / r * Kokkos::sin(theta);
            } else {
                //Compute the (2,2) coefficient of the inverse Jacobian matrix.
                return 1 / r * Kokkos::cos(theta);
            }
        }
    }


    CartesianToCircular<X, Y, R, Theta> get_inverse_mapping() const
    {
        return CartesianToCircular<X, Y, R, Theta>();
    }
};

namespace mapping_detail {
template <class X, class Y, class R, class Theta, class ExecSpace>
struct MappingAccessibility<ExecSpace, CircularToCartesian<R, Theta, X, Y>> : std::true_type
{
};

template <class X, class Y, class R, class Theta>
struct IsCurvilinear2DMapping<CircularToCartesian<R, Theta, X, Y>> : std::true_type
{
};

template <class X, class Y, class R, class Theta>
struct SingularOPointInvJacobian<CircularToCartesian<R, Theta, X, Y>> : std::true_type
{
};

} // namespace mapping_detail
```


