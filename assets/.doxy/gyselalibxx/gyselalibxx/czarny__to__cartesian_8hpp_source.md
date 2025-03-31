

# File czarny\_to\_cartesian.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**czarny\_to\_cartesian.hpp**](czarny__to__cartesian_8hpp.md)

[Go to the documentation of this file](czarny__to__cartesian_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "tensor.hpp"
#include "view.hpp"

// Pre-declaration of analytical inverse
template <class X, class Y, class R, class Theta>
class CartesianToCzarny;

template <class R, class Theta, class X, class Y>
class CzarnyToCartesian
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

private:
    double m_epsilon;
    double m_e;

public:
    CzarnyToCartesian(double epsilon, double e) : m_epsilon(epsilon), m_e(e) {}

    KOKKOS_FUNCTION CzarnyToCartesian(CzarnyToCartesian const& other)
        : m_epsilon(other.epsilon())
        , m_e(other.e())
    {
    }

    CzarnyToCartesian(CzarnyToCartesian&& x) = default;

    ~CzarnyToCartesian() = default;

    CzarnyToCartesian& operator=(CzarnyToCartesian const& x) = default;

    CzarnyToCartesian& operator=(CzarnyToCartesian&& x) = default;

    KOKKOS_FUNCTION double epsilon() const
    {
        return m_epsilon;
    }

    KOKKOS_FUNCTION double e() const
    {
        return m_e;
    }

    KOKKOS_FUNCTION Coord<X, Y> operator()(Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        const double tmp1
                = Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * Kokkos::cos(theta)) + 1.0);

        const double x = (1.0 - tmp1) / m_epsilon;
        const double y = m_e * r * Kokkos::sin(theta)
                         / (Kokkos::sqrt(1.0 - 0.25 * m_epsilon * m_epsilon) * (2.0 - tmp1));

        return Coord<X, Y>(x, y);
    }

    KOKKOS_FUNCTION double jacobian(Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        const double xi = Kokkos::sqrt(1. / (1. - m_epsilon * m_epsilon * 0.25));
        return -r / Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * Kokkos::cos(theta))) * m_e
               * xi
               / (2 - Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * Kokkos::cos(theta))));
    }

    KOKKOS_FUNCTION DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix(
            Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);

        const double sin_theta = Kokkos::sin(theta);
        const double cos_theta = Kokkos::cos(theta);
        const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
        const double xi = Kokkos::sqrt(xi2);
        const double sqrt_eps = Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * cos_theta) + 1.0);
        const double sqrt_eps_2 = 2.0 - sqrt_eps;

        DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix;
        ddcHelper::get<X, R_cov>(jacobian_matrix) = -cos_theta / sqrt_eps;
        ddcHelper::get<X, Theta_cov>(jacobian_matrix) = r * sin_theta / sqrt_eps;
        ddcHelper::get<Y, R_cov>(jacobian_matrix) = m_e * m_epsilon * r * sin_theta * cos_theta * xi
                                                            / (sqrt_eps_2 * sqrt_eps_2 * sqrt_eps)
                                                    + m_e * sin_theta * xi / sqrt_eps_2;
        ddcHelper::get<Y, Theta_cov>(jacobian_matrix)
                = r
                  * (-m_e * m_epsilon * r * sin_theta * sin_theta * xi
                             / (sqrt_eps_2 * sqrt_eps_2 * sqrt_eps)
                     + m_e * cos_theta * xi / sqrt_eps_2);
        return jacobian_matrix;
    }



    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double jacobian_component(Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (1,1), i.e dx/dr
            return -Kokkos::cos(theta)
                   / Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * Kokkos::cos(theta)) + 1.0);
        } else if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, Theta_cov>) {
            // Component (1,2), i.e dx/dtheta
            return r * Kokkos::sin(theta)
                   / Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * Kokkos::cos(theta)) + 1.0);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (2,1), i.e dy/dr
            const double sin_theta = Kokkos::sin(theta);
            const double cos_theta = Kokkos::cos(theta);
            const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
            const double xi = Kokkos::sqrt(xi2);
            const double tmp1 = Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * cos_theta) + 1.0);
            const double tmp2 = 2.0 - tmp1;
            return m_e * m_epsilon * r * sin_theta * cos_theta * xi / (tmp2 * tmp2 * tmp1)
                   + m_e * sin_theta * xi / tmp2;
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, Theta_cov>) {
            // Component (2,2), i.e dy/dtheta
            const double sin_theta = Kokkos::sin(theta);
            const double cos_theta = Kokkos::cos(theta);
            const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
            const double xi = Kokkos::sqrt(xi2);
            const double tmp1 = Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * cos_theta) + 1.0);
            const double tmp2 = 2.0 - tmp1;
            return r
                   * (-m_e * m_epsilon * r * sin_theta * sin_theta * xi / (tmp2 * tmp2 * tmp1)
                      + m_e * cos_theta * xi / tmp2);
        }
    }


    KOKKOS_FUNCTION DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X_cov, Y_cov>>
    inv_jacobian_matrix(Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);

        assert(r >= 1e-15);

        const double sin_theta = Kokkos::sin(theta);
        const double cos_theta = Kokkos::cos(theta);
        const double xi = Kokkos::sqrt(1. / (1. - m_epsilon * m_epsilon * 0.25));
        const double divisor = 2 - Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * cos_theta));

        const double fact_1 = 1 / Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * cos_theta));
        const double fact_2 = m_e * m_epsilon * xi * r * sin_theta * fact_1 / divisor / divisor;
        const double fact_3 = m_e * xi / divisor;

        DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X_cov, Y_cov>> matrix;
        ddcHelper::get<R, X_cov>(matrix)
                = -1 / fact_1 * (-sin_theta * fact_2 + cos_theta * fact_3) / fact_3;
        ddcHelper::get<R, Y_cov>(matrix) = sin_theta / fact_3;
        ddcHelper::get<Theta, X_cov>(matrix)
                = 1 / r / fact_1 * (cos_theta * fact_2 + sin_theta * fact_3) / fact_3;
        ddcHelper::get<Theta, Y_cov>(matrix) = 1 / r * cos_theta / fact_3;
        return matrix;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double inv_jacobian_component(Coord<R, Theta> const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<R, Theta>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<X_cov, Y_cov>>);

        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        const double cos_theta = Kokkos::cos(theta);
        const double xi = Kokkos::sqrt(1. / (1. - m_epsilon * m_epsilon * 0.25));
        const double divisor = 2 - Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * cos_theta));

        if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, X_cov>) {
            //Compute the (1,1) coefficient of the inverse Jacobian matrix
            const double sin_theta = Kokkos::sin(theta);

            const double fact_1
                    = 1 / Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * cos_theta));
            const double fact_2 = m_e * m_epsilon * xi * r * sin_theta * fact_1 / divisor / divisor;
            const double fact_3 = m_e * xi / divisor;

            return -1 / fact_1 * (-sin_theta * fact_2 + cos_theta * fact_3) / fact_3;
        } else if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, Y_cov>) {
            //Compute the (1,2) coefficient of the Jacobian matrix, i.e J^x_theta.
            const double sin_theta = Kokkos::sin(theta);

            const double fact_3 = m_e * xi / divisor;
            return sin_theta / fact_3;
        } else if constexpr (std::is_same_v<IndexTag1, Theta> && std::is_same_v<IndexTag2, X_cov>) {
            //Compute the (2,1) coefficient of the inverse Jacobian matrix.
            assert(r >= 1e-15);

            const double sin_theta = Kokkos::sin(theta);
            const double fact_1
                    = 1 / Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * cos_theta));
            const double fact_2 = m_e * m_epsilon * xi * r * sin_theta * fact_1 / divisor / divisor;
            const double fact_3 = m_e * xi / divisor;

            return 1 / r / fact_1 * (cos_theta * fact_2 + sin_theta * fact_3) / fact_3;
        } else {
            //Compute the (2,2) coefficient of the inverse Jacobian matrix.
            assert(r >= 1e-15);

            const double fact_3 = m_e * xi / divisor;
            return 1 / r * cos_theta / fact_3;
        }
    }


    CartesianToCzarny<X, Y, R, Theta> get_inverse_mapping() const
    {
        return CartesianToCzarny<X, Y, R, Theta>(m_epsilon, m_e);
    }
};

namespace mapping_detail {
template <class X, class Y, class R, class Theta, class ExecSpace>
struct MappingAccessibility<ExecSpace, CzarnyToCartesian<R, Theta, X, Y>> : std::true_type
{
};

template <class X, class Y, class R, class Theta>
struct IsCurvilinear2DMapping<CzarnyToCartesian<X, Y, R, Theta>> : std::true_type
{
};

template <class X, class Y, class R, class Theta>
struct SingularOPointInvJacobian<CzarnyToCartesian<R, Theta, X, Y>> : std::true_type
{
};

} // namespace mapping_detail
```


