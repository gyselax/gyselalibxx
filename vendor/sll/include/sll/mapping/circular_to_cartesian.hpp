#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include "coordinate_converter.hpp"
#include "curvilinear2d_to_cartesian.hpp"
#include "jacobian.hpp"
#include "pseudo_cartesian_compatible_mapping.hpp"

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
template <class X, class Y, class R, class Theta>
class CircularToCartesian
    : public CoordinateConverter<ddc::Coordinate<X, Y>, ddc::Coordinate<R, Theta>>
    , public CoordinateConverter<ddc::Coordinate<R, Theta>, ddc::Coordinate<X, Y>>
    , public Jacobian<ddc::Coordinate<R, Theta>>
    , public PseudoCartesianCompatibleMapping
    , public Curvilinear2DToCartesian<X, Y, R, Theta>
{
public:
    /// @brief Indicate the first physical coordinate.
    using cartesian_tag_x = typename Curvilinear2DToCartesian<X, Y, R, Theta>::cartesian_tag_x;
    /// @brief Indicate the second physical coordinate.
    using cartesian_tag_y = typename Curvilinear2DToCartesian<X, Y, R, Theta>::cartesian_tag_y;
    /// @brief Indicate the first logical coordinate.
    using curvilinear_tag_r = typename Curvilinear2DToCartesian<X, Y, R, Theta>::curvilinear_tag_r;
    /// @brief Indicate the second logical coordinate.
    using curvilinear_tag_theta =
            typename Curvilinear2DToCartesian<X, Y, R, Theta>::curvilinear_tag_theta;
    /// The type of the Jacobian matrix and its inverse
    using Matrix_2x2 = typename Jacobian<ddc::Coordinate<R, Theta>>::Matrix_2x2;

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
     * @brief Instantiate a Curvilinear2DToCartesian from another temporary CircularToCartesian (rvalue).
     *
     * @param[in] x
     * 		Curvilinear2DToCartesian mapping used to instantiate the new one.
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

    KOKKOS_FUNCTION ddc::Coordinate<X, Y> operator()(
            ddc::Coordinate<R, Theta> const& coord) const final
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        const double x = r * Kokkos::cos(theta);
        const double y = r * Kokkos::sin(theta);
        return ddc::Coordinate<X, Y>(x, y);
    }

    KOKKOS_FUNCTION ddc::Coordinate<R, Theta> operator()(
            ddc::Coordinate<X, Y> const& coord) const final
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        const double r = Kokkos::sqrt(x * x + y * y);
        const double theta = Kokkos::atan2(y, x);
        return ddc::Coordinate<R, Theta>(r, theta);
    }

    double jacobian(ddc::Coordinate<R, Theta> const& coord) const final
    {
        double r = ddc::get<R>(coord);
        return r;
    }


    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given indendently with the functions
     * jacobian_11, jacobian_12,  jacobian_21 and jacobian_22.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @param[out] matrix
     * 				The Jacobian matrix returned.
     */
    void jacobian_matrix(ddc::Coordinate<R, Theta> const& coord, Matrix_2x2& matrix) const final
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        matrix[0][0] = std::cos(theta);
        matrix[0][1] = -r * std::sin(theta);
        matrix[1][0] = std::sin(theta);
        matrix[1][1] = r * std::cos(theta);
    }

    double jacobian_11(ddc::Coordinate<R, Theta> const& coord) const final
    {
        const double theta = ddc::get<Theta>(coord);
        return std::cos(theta);
    }

    double jacobian_12(ddc::Coordinate<R, Theta> const& coord) const final
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        return -r * std::sin(theta);
    }

    double jacobian_21(ddc::Coordinate<R, Theta> const& coord) const final
    {
        const double theta = ddc::get<Theta>(coord);
        return std::sin(theta);
    }

    double jacobian_22(ddc::Coordinate<R, Theta> const& coord) const final
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        return r * std::cos(theta);
    }


    /**
     * @brief Compute full inverse Jacobian matrix.
     *
     * For some computations, we need the complete inverse Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given indendently with the functions
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
    void inv_jacobian_matrix(ddc::Coordinate<R, Theta> const& coord, Matrix_2x2& matrix) const final
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        assert(fabs(r) >= 1e-15);
        matrix[0][0] = std::cos(theta);
        matrix[0][1] = std::sin(theta);
        matrix[1][0] = -1 / r * std::sin(theta);
        matrix[1][1] = 1 / r * std::cos(theta);
    }

    double inv_jacobian_11(ddc::Coordinate<R, Theta> const& coord) const final
    {
        const double theta = ddc::get<Theta>(coord);
        return std::cos(theta);
    }

    double inv_jacobian_12(ddc::Coordinate<R, Theta> const& coord) const final
    {
        const double theta = ddc::get<Theta>(coord);
        return std::sin(theta);
    }

    double inv_jacobian_21(ddc::Coordinate<R, Theta> const& coord) const final
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        assert(fabs(r) >= 1e-15);
        return -1 / r * std::sin(theta);
    }

    double inv_jacobian_22(ddc::Coordinate<R, Theta> const& coord) const final
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        assert(fabs(r) >= 1e-15);
        return 1 / r * std::cos(theta);
    }



    /**
     * @brief  Compute the full Jacobian matrix from the mapping to the pseudo-Cartesian mapping at the central point.
     *
     *
     * Here, as @f$ \mathcal{G} =  \mathcal{F} @f$ (see PseudoCartesianCompatibleMapping), the Jacobian matrix of
     * @f$(\mathcal{F} \circ \mathcal{G}^{-1})^{-1} @f$ is the identity matrix.
     * So, the pseudo-Cartesian Jacobian matrix for a circular mapping is given by :
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{11}(0, \theta) = 1, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{12}(0, \theta) = 0, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{21}(0, \theta) = 0, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{22}(0, \theta) = 1. @f$
     *
     *
     * @param[out] matrix
     *      The pseudo-Cartesian matrix evaluated at the central point.
     *
     *
     * @see DiscreteToCartesian
     * @see BslAdvection
     * @see AdvectionDomain
     */
    void to_pseudo_cartesian_jacobian_center_matrix(Matrix_2x2& matrix) const final
    {
        matrix[0][0] = 1.;
        matrix[0][1] = 0.;
        matrix[1][0] = 0.;
        matrix[1][1] = 1.;
    }

    /**
     * @brief Compute the (1,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{11}(0, \theta) = 1. @f$
     *
     * @return A double with the (1,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see to_pseudo_cartesian_jacobian_center_matrix
     */
    double to_pseudo_cartesian_jacobian_11_center() const final
    {
        return 1.;
    }

    /**
     * @brief Compute the (1,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{12}(0, \theta) = 0. @f$
     *
     * @return A double with the (1,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see to_pseudo_cartesian_jacobian_center_matrix
     */
    double to_pseudo_cartesian_jacobian_12_center() const final
    {
        return 0.;
    }

    /**
     * @brief Compute the (2,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{21}(0, \theta) = 0. @f$
     *
     * @return A double with the (2,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see to_pseudo_cartesian_jacobian_center_matrix
     */
    double to_pseudo_cartesian_jacobian_21_center() const final
    {
        return 0.;
    }

    /**
     * @brief Compute the (2,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{22}(0, \theta) = 1. @f$
     *
     * @return A double with the (2,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see to_pseudo_cartesian_jacobian_center_matrix
     */
    double to_pseudo_cartesian_jacobian_22_center() const final
    {
        return 1.;
    }
};
