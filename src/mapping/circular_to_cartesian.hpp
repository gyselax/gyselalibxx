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

/**
 * @brief A class for describing the circular 2D mapping.
 *
 * The mapping @f$ (r,\theta)\mapsto (x,y) @f$ is defined as follow :
 *
 * @f$ x(r,\theta) = r \cos(\theta) + x_0,@f$
 *
 * @f$ y(r,\theta) = r \sin(\theta) + y_0.@f$
 *
 * It and its Jacobian matrix are invertible everywhere except for @f$ r = 0 @f$.
 *
 * The Jacobian matrix coefficients are defined as follow
 *
 * @f$ J_{11}(r,\theta)  = \cos(\theta)@f$
 *
 * @f$ J_{12}(r,\theta)  = - r \sin(\theta)@f$
 *
 * @f$ J_{21}(r,\theta)  = \sin(\theta)@f$
 *
 * @f$ J_{22}(r,\theta)  = r \cos(\theta)@f$
 *
 * and the matrix determinant: @f$ det(J) = r @f$.
 *
 */
template <class R, class Theta, class X, class Y>
class CircularToCartesian
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
    using CoordArg = Coord<R, Theta>;
    /// The type of the result of the function described by this mapping
    using CoordResult = Coord<X, Y>;

    /// @brief The covariant form of the first physical coordinate.
    using X_cov = typename X::Dual;
    /// @brief The covariant form of the second physical coordinate.
    using Y_cov = typename Y::Dual;
    /// @brief The covariant form of the first logical coordinate.
    using R_cov = typename R::Dual;
    /// @brief The covariant form of the second logical coordinate.
    using Theta_cov = typename Theta::Dual;

private:
    double m_x0;
    double m_y0;

public:
    /**
     * @brief Instantiate a CircularToCartesian from parameters.
     *
     * @param[in] x0 The x-coordinate of the centre of the circle (0 by default).
     * @param[in] y0 The y-coordinate of the centre of the circle (0 by default).
     */
    explicit KOKKOS_FUNCTION CircularToCartesian(double x0 = 0.0, double y0 = 0.0)
        : m_x0(x0)
        , m_y0(y0)
    {
    }

    /**
     * @brief Instantiate a CircularToCartesian from another CircularToCartesian (lvalue).
     *
     * @param[in] other
     * 		CircularToCartesian mapping used to instantiate the new one.
     */
    KOKKOS_DEFAULTED_FUNCTION CircularToCartesian(CircularToCartesian const& other) = default;

    /**
     * @brief Instantiate a CircularToCartesian from another temporary CircularToCartesian (rvalue).
     *
     * @param[in] x
     * 		CircularToCartesian mapping used to instantiate the new one.
     */
    CircularToCartesian(CircularToCartesian&& x) = default;

    KOKKOS_DEFAULTED_FUNCTION ~CircularToCartesian() = default;

    /**
     * @brief Assign a CircularToCartesian from another CircularToCartesian (lvalue).
     *
     * @param[in] x
     * 		CircularToCartesian mapping used to assign.
     *
     * @return The CircularToCartesian assigned.
     */
    CircularToCartesian& operator=(CircularToCartesian const& x) = default;

    /**
     * @brief Assign a CircularToCartesian from another temporary CircularToCartesian (rvalue).
     *
     * @param[in] x
     * 		CircularToCartesian mapping used to assign.
     *
     * @return The CircularToCartesian assigned.
     */
    CircularToCartesian& operator=(CircularToCartesian&& x) = default;

    /**
     * @brief Convert the @f$ (r, \theta) @f$ coordinate to the equivalent (x,y) coordinate.
     *
     * @param[in] coord The coordinate to be converted.
     *
     * @return The equivalent coordinate.
     */
    KOKKOS_FUNCTION Coord<X, Y> operator()(Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        const double x = r * Kokkos::cos(theta) + m_x0;
        const double y = r * Kokkos::sin(theta) + m_y0;
        return Coord<X, Y>(x, y);
    }

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian(Coord<R, Theta> const& coord) const
    {
        double r = ddc::get<R>(coord);
        return r;
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

    /**
     * @brief Compute the (i,j) coefficient of the Jacobian matrix.
     * 
     * @param[in] coord
     *              The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (i,j) coefficient of the Jacobian matrix.
     */
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


    /**
     * @brief Compute full inverse Jacobian matrix.
     *
     * For some computations, we need the complete inverse Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given independently with the function inv_jacobian_component.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The inverse Jacobian matrix.
     */
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


    /**
     * @brief Compute the (i,j) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the centre point.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (i,j) coefficient of the inverse Jacobian matrix.
     */
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


    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
    KOKKOS_INLINE_FUNCTION CartesianToCircular<X, Y, R, Theta> get_inverse_mapping() const
    {
        return CartesianToCircular<X, Y, R, Theta>(m_x0, m_y0);
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
