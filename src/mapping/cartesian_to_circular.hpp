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

/**
 * @brief A class for describing the circular 2D mapping.
 *
 * The mapping @f$ (x,y)\mapsto (r,\theta) @f$ is defined as follow :
 *
 * @f$ r(x,y) = \sqrt x^2+y^2 ,@f$
 *
 * @f$ \theta(x,y) = atan2(\frac{y}{x}) .@f$
 *
 * It and its Jacobian matrix are invertible everywhere except for @f$ r = 0 @f$.
 *
 * The Jacobian matrix coefficients are defined as follow
 *
 * @f$ J_{11}(r,\theta)  =\frac{2x}{\sqrt{x^2+y^2}}  @f$
 *
 * @f$ J_{12}(r,\theta)  =\frac{2y}{\sqrt{x^2+y^2}}  @f$
 *
 * @f$ J_{21}(r,\theta)  =\frac{-y}{(x^2+y^2)^2}  @f$
 *
 * @f$ J_{22}(r,\theta)  =\frac{x}{(x^2+y^2)^2}  @f$
 *
 * and the matrix determinant: @f$ det(J) = r @f$.
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

    /// @brief The covariant form of the first physical coordinate.
    using X_cov = typename X::Dual;
    /// @brief The covariant form of the second physical coordinate.
    using Y_cov = typename Y::Dual;
    /// @brief The covariant form of the first logical coordinate.
    using R_cov = typename R::Dual;
    /// @brief The covariant form of the second logical coordinate.
    using Theta_cov = typename Theta::Dual;

public:
    CartesianToCircular() = default;

    /**
     * @brief Instantiate a CartesianToCircular from another CartesianToCircular (lvalue).
     *
     * @param[in] other
     * 		CartesianToCircular mapping used to instantiate the new one.
     */
    KOKKOS_FUNCTION CartesianToCircular(CartesianToCircular const& other) {}

    /**
     * @brief Instantiate a Curvilinear2DToCartesian from another temporary CartesianToCircular (rvalue).
     *
     * @param[in] x
     * 		Curvilinear2DToCartesian mapping used to instantiate the new one.
     */
    CartesianToCircular(CartesianToCircular&& x) = default;

    ~CartesianToCircular() = default;

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
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        const double r = Kokkos::sqrt(x * x + y * y);
        const double theta = Kokkos::atan2(y, x);
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
    KOKKOS_FUNCTION double jacobian(Coord<X, Y> const& coord)
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        return 2. * (x * x - y * y) / Kokkos::pow(x * x + y * y, 1.5);
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

    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
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
