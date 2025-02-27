// SPDX-License-Identifier: MIT
#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "view.hpp"

// Pre-declaration of analytical inverse
template <class X, class Y, class R, class Theta>
class CartesianToCircular;

/**
 * @brief A class for describing the circular 2D mapping.
 *
 * The mapping @f$ (r,\theta)\mapsto (x,y) @f$ is defined as follow :
 *
 * @f$ x(r,\theta) = r \cos(\theta),@f$
 *
 * @f$ y(r,\theta) = r \sin(\theta).@f$
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

public:
    CircularToCartesian() = default;

    /**
     * @brief Instantiate a CircularToCartesian from another CircularToCartesian (lvalue).
     *
     * @param[in] other
     * 		CircularToCartesian mapping used to instantiate the new one.
     */
    KOKKOS_FUNCTION CircularToCartesian(CircularToCartesian const& other) {}

    /**
     * @brief Instantiate a CircularToCartesian from another temporary CircularToCartesian (rvalue).
     *
     * @param[in] x
     * 		CircularToCartesian mapping used to instantiate the new one.
     */
    CircularToCartesian(CircularToCartesian&& x) = default;

    ~CircularToCartesian() = default;

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
        const double x = r * Kokkos::cos(theta);
        const double y = r * Kokkos::sin(theta);
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
     * The coefficients can be given independently with the functions
     * jacobian_11, jacobian_12,  jacobian_21 and jacobian_22.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @param[out] matrix
     * 				The Jacobian matrix returned.
     */
    KOKKOS_FUNCTION void jacobian_matrix(Coord<R, Theta> const& coord, Matrix_2x2& matrix) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        matrix[0][0] = Kokkos::cos(theta);
        matrix[0][1] = -r * Kokkos::sin(theta);
        matrix[1][0] = Kokkos::sin(theta);
        matrix[1][1] = r * Kokkos::cos(theta);
    }

    /**
     * @brief Compute the (1,1) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (1,1) coefficient of the Jacobian matrix is given by @f$ \frac{\partial x}{\partial r} @f$.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (1,1) coefficient of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian_11(Coord<R, Theta> const& coord) const
    {
        const double theta = ddc::get<Theta>(coord);
        return Kokkos::cos(theta);
    }

    /**
     * @brief Compute the (1,2) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (1,2) coefficient of the Jacobian matrix is given by @f$ \frac{\partial x}{\partial \theta} @f$.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (1,2) coefficient of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian_12(Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        return -r * Kokkos::sin(theta);
    }

    /**
     * @brief Compute the (2,1) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (2,1) coefficient of the Jacobian matrix is given by @f$ \frac{\partial y}{\partial r} @f$.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the Jacobian matrix. .
     *
     * @return A double with the value of the (2,1) coefficient of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian_21(Coord<R, Theta> const& coord) const
    {
        const double theta = ddc::get<Theta>(coord);
        return Kokkos::sin(theta);
    }

    /**
     * @brief Compute the (2,2) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (2,2) coefficient of the Jacobian matrix is given by @f$ \frac{\partial y}{\partial \theta} @f$.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (2,2) coefficient of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian_22(Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        return r * Kokkos::cos(theta);
    }


    /**
     * @brief Compute full inverse Jacobian matrix.
     *
     * For some computations, we need the complete inverse Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given independently with the functions
     * inv_jacobian_11, inv_jacobian_12, inv_jacobian_21 and inv_jacobian_22.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @param[out] matrix
     * 				The inverse Jacobian matrix returned.
     *
     *
     * @see Jacobian::inv_jacobian_11
     * @see Jacobian::inv_jacobian_12
     * @see Jacobian::inv_jacobian_21
     * @see Jacobian::inv_jacobian_22
     */
    KOKKOS_FUNCTION void inv_jacobian_matrix(Coord<R, Theta> const& coord, Matrix_2x2& matrix) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        assert(fabs(r) >= 1e-15);
        matrix[0][0] = Kokkos::cos(theta);
        matrix[0][1] = Kokkos::sin(theta);
        matrix[1][0] = -1 / r * Kokkos::sin(theta);
        matrix[1][1] = 1 / r * Kokkos::cos(theta);
    }

    /**
     * @brief Compute the (1,1) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the centre point.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (1,1) coefficient of the inverse Jacobian matrix.
     */
    KOKKOS_FUNCTION double inv_jacobian_11(Coord<R, Theta> const& coord) const
    {
        const double theta = ddc::get<Theta>(coord);
        return Kokkos::cos(theta);
    }

    /**
     * @brief Compute the (1,2) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the centre point.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (1,2) coefficient of the inverse Jacobian matrix.
     */
    KOKKOS_FUNCTION double inv_jacobian_12(Coord<R, Theta> const& coord) const
    {
        const double theta = ddc::get<Theta>(coord);
        return Kokkos::sin(theta);
    }

    /**
     * @brief Compute the (2,1) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the centre point.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (2,1) coefficient of the inverse Jacobian matrix.
     */
    KOKKOS_FUNCTION double inv_jacobian_21(Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        assert(fabs(r) >= 1e-15);
        return -1 / r * Kokkos::sin(theta);
    }

    /**
     * @brief Compute the (2,2) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the centre point.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (2,2) coefficient of the inverse Jacobian matrix.
     */
    KOKKOS_FUNCTION double inv_jacobian_22(Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        assert(fabs(r) >= 1e-15);
        return 1 / r * Kokkos::cos(theta);
    }

    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
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
