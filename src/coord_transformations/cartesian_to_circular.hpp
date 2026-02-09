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

/**
 * @brief A class for describing the circular 2D mapping.
 *
 * The mapping @f$ (x,y)\mapsto (r,\theta) @f$ is defined as follow :
 *
 * @f$ r(x,y) = \sqrt (x-x_0)^2+(y-y_0)^2 ,@f$
 *
 * @f$ \theta(x,y) = atan2(\frac{y-y_0}{x-x_0}) .@f$
 *
 * It and its Jacobian matrix are invertible everywhere except for @f$ r = 0 @f$.
 *
 * The Jacobian matrix coefficients are defined as follow
 *
 * @f$ J_{11}(x,y)  =\frac{x-x_0}{\sqrt{(x-x_0)^2+(y-y_0)^2}}  @f$
 *
 * @f$ J_{12}(x,y)  =\frac{y-y_0}{\sqrt{(x-x_0)^2+(y-y_0)^2}}  @f$
 *
 * @f$ J_{21}(x,y)  =\frac{-(y-y_0)}{(x-x_0)^2+(y-y_0)^2}  @f$
 *
 * @f$ J_{22}(x,y)  =\frac{x-x_0}{(x-x_0)^2+(y-y_0)^2}  @f$
 *
 * and the matrix determinant: @f$ det(J) = 1/((x-x_0)^2+(y-y_0)^2) @f$.
 *
 */
template <class X, class Y, class R, class Theta>
class CartesianToCircular
{
public:
    /// @brief Indicate the first physical coordinate.
    using cartesian_tag_x = X;
    /// @brief Indicate the second physical coordinate.
    using cartesian_tag_y = Y;
    /// @brief Indicate the first logical coordinate.
    using curvilinear_tag_r = R;
    /// @brief Indicate the second logical coordinate.
    using curvilinear_tag_theta = Theta;

    /// The type of the argument of the function described by this mapping
    using CoordArg = Coord<X, Y>;
    /// The type of the result of the function described by this mapping
    using CoordResult = Coord<R, Theta>;
    /// The type of the coordinate that can be used to evaluate the Jacobian of this mapping
    using CoordJacobian = CoordArg;

    /// @brief The covariant form of the first physical coordinate.
    using X_cov = typename X::Dual;
    /// @brief The covariant form of the second physical coordinate.
    using Y_cov = typename Y::Dual;
    /// @brief The covariant form of the first logical coordinate.
    using R_cov = typename R::Dual;
    /// @brief The covariant form of the second logical coordinate.
    using Theta_cov = typename Theta::Dual;

private:
    Coord<X, Y> m_o_point;

public:
    /**
     * @brief Instantiate a CartesianToCircular from parameters.
     *
     * @param[in] o_point The (x,y)-coordinate of the centre of the circle ((0,0) by default).
     */
    explicit KOKKOS_FUNCTION CartesianToCircular(Coord<X, Y> o_point = Coord<X, Y>(0.0, 0.0))
        : m_o_point(o_point)
    {
    }

    /**
     * @brief Instantiate a CartesianToCircular from another CartesianToCircular (lvalue).
     *
     * @param[in] other
     * 		CartesianToCircular mapping used to instantiate the new one.
     */
    KOKKOS_DEFAULTED_FUNCTION CartesianToCircular(CartesianToCircular const& other) = default;

    /**
     * @brief Instantiate a Curvilinear2DToCartesian from another temporary CartesianToCircular (rvalue).
     *
     * @param[in] x
     * 		Curvilinear2DToCartesian mapping used to instantiate the new one.
     */
    CartesianToCircular(CartesianToCircular&& x) = default;

    KOKKOS_DEFAULTED_FUNCTION ~CartesianToCircular() = default;

    /**
     * @brief Assign a CartesianToCircular from another CartesianToCircular (lvalue).
     *
     * @param[in] x
     * 		CartesianToCircular mapping used to assign.
     *
     * @return The CartesianToCircular assigned.
     */
    CartesianToCircular& operator=(CartesianToCircular const& x) = default;

    /**
     * @brief Assign a CartesianToCircular from another temporary CartesianToCircular (rvalue).
     *
     * @param[in] x
     * 		CartesianToCircular mapping used to assign.
     *
     * @return The CartesianToCircular assigned.
     */
    CartesianToCircular& operator=(CartesianToCircular&& x) = default;

    /**
     * @brief Convert the coordinate (x,y) to the equivalent @f$ (r, \theta) @f$ coordinate.
     *
     * @param[in] coord The coordinate to be converted.
     *
     * @return The equivalent coordinate.
     */
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

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian(Coord<X, Y> const& coord) const
    {
        const double x = ddc::get<X>(coord) - ddc::get<X>(m_o_point);
        const double y = ddc::get<Y>(coord) - ddc::get<Y>(m_o_point);
        return 1. / Kokkos::sqrt(x * x + y * y);
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given independently with the function jacobian_component.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The Jacobian matrix.
     */
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

    /**
     * @brief Compute the (i,j) coefficient of the Jacobian matrix.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (1,1) coefficient of the Jacobian matrix.
     */
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

    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
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
