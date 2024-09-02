#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include "analytical_invertible_curvilinear2d_to_cartesian.hpp"

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
 *
 * @see AnalyticalInvertibleCurvilinear2DToCartesian
 */
template <class X, class Y, class R, class Theta>
class CircularToCartesian : public AnalyticalInvertibleCurvilinear2DToCartesian<X, Y, R, Theta>
{
public:
    /**
     * @brief Indicate the first physical coordinate.
     */
    using cartesian_tag_x = X;
    /**
     * @brief Indicate the second physical coordinate.
     */
    using cartesian_tag_y = Y;
    /**
     * @brief Indicate the first logical coordinate.
     */
    using circular_tag_r = R;
    /**
     * @brief Indicate the second logical coordinate.
     */
    using circular_tag_theta = Theta;

    /**
     * @brief Define a 2x2 matrix with an 2D array of an 2D array.
     */
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

public:
    CircularToCartesian() = default;

    /**
     * @brief Instantiate a CircularToCartesian from another CircularToCartesian (lvalue).
     *
     * @param[in] other
     * 		CircularToCartesian mapping used to instantiate the new one.
     */
    CircularToCartesian(CircularToCartesian const& other) = default;

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

    ddc::Coordinate<X, Y> operator()(ddc::Coordinate<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        const double x = r * std::cos(theta);
        const double y = r * std::sin(theta);
        return ddc::Coordinate<X, Y>(x, y);
    }

    ddc::Coordinate<R, Theta> operator()(ddc::Coordinate<X, Y> const& coord) const
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        const double r = std::sqrt(x * x + y * y);
        const double theta = std::atan2(y, x);
        return ddc::Coordinate<R, Theta>(r, theta);
    }

    double jacobian(ddc::Coordinate<R, Theta> const& coord) const final
    {
        double r = ddc::get<R>(coord);
        return r;
    }


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
     * Here, as @f$ \mathcal{G} =  \mathcal{F} @f$ (see DiscreteToCartesian), the Jacobian matrix of
     * @f$(\mathcal{F} \circ \mathcal{G}^{-1})^{-1} @f$ is the identity matrix.
     * So, the pseudo-Cartesian Jacobian matrix for a circular mapping is given by :
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{11}(0, \theta) = 1, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{12}(0, \theta) = 0, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{21}(0, \theta) = 0, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{22}(0, \theta) = 1. @f$
     *
     *
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     * @param[out] matrix
     *      The pseudo-Cartesian matrix evaluated at the central point.
     *
     *
     * @see DiscreteToCartesian
     * @see BslAdvection
     * @see AdvectionDomain
     */
    template <class IdxRange>
    void to_pseudo_cartesian_jacobian_center_matrix(IdxRange const& grid, Matrix_2x2& matrix) const
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
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (1,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see CircularToCartesian::to_pseudo_cartesian_jacobian_center_matrix
     */
    template <class IdxRange>
    double to_pseudo_cartesian_jacobian_11_center(IdxRange const& grid) const
    {
        return 1.;
    }

    /**
     * @brief Compute the (1,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{12}(0, \theta) = 0. @f$
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (1,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see CircularToCartesian::to_pseudo_cartesian_jacobian_center_matrix
     */
    template <class IdxRange>
    double to_pseudo_cartesian_jacobian_12_center(IdxRange const& grid) const
    {
        return 0.;
    }

    /**
     * @brief Compute the (2,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{21}(0, \theta) = 0. @f$
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (2,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see CircularToCartesian::to_pseudo_cartesian_jacobian_center_matrix
     */
    template <class IdxRange>
    double to_pseudo_cartesian_jacobian_21_center(IdxRange const& grid) const
    {
        return 0.;
    }

    /**
     * @brief Compute the (2,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{22}(0, \theta) = 1. @f$
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (2,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see CircularToCartesian::to_pseudo_cartesian_jacobian_center_matrix
     */
    template <class IdxRange>
    double to_pseudo_cartesian_jacobian_22_center(IdxRange const& grid) const
    {
        return 1.;
    }
};
