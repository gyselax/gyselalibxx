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
template <class DimX, class DimY, class DimR, class DimP>
class CircularToCartesian
    : public AnalyticalInvertibleCurvilinear2DToCartesian<DimX, DimY, DimR, DimP>
{
public:
    /**
     * @brief Indicate the first physical coordinate.
     */
    using cartesian_tag_x = DimX;
    /**
     * @brief Indicate the second physical coordinate.
     */
    using cartesian_tag_y = DimY;
    /**
     * @brief Indicate the first logical coordinate.
     */
    using circular_tag_r = DimR;
    /**
     * @brief Indicate the second logical coordinate.
     */
    using circular_tag_p = DimP;

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

    ddc::Coordinate<DimX, DimY> operator()(ddc::Coordinate<DimR, DimP> const& coord) const
    {
        const double r = ddc::get<DimR>(coord);
        const double p = ddc::get<DimP>(coord);
        const double x = r * std::cos(p);
        const double y = r * std::sin(p);
        return ddc::Coordinate<DimX, DimY>(x, y);
    }

    ddc::Coordinate<DimR, DimP> operator()(ddc::Coordinate<DimX, DimY> const& coord) const
    {
        const double x = ddc::get<DimX>(coord);
        const double y = ddc::get<DimY>(coord);
        const double r = std::sqrt(x * x + y * y);
        const double p = std::atan2(y, x);
        return ddc::Coordinate<DimR, DimP>(r, p);
    }

    double jacobian(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        double r = ddc::get<DimR>(coord);
        return r;
    }


    void jacobian_matrix(ddc::Coordinate<DimR, DimP> const& coord, Matrix_2x2& matrix) const final
    {
        const double r = ddc::get<DimR>(coord);
        const double p = ddc::get<DimP>(coord);
        matrix[0][0] = std::cos(p);
        matrix[0][1] = -r * std::sin(p);
        matrix[1][0] = std::sin(p);
        matrix[1][1] = r * std::cos(p);
    }

    double jacobian_11(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double p = ddc::get<DimP>(coord);
        return std::cos(p);
    }

    double jacobian_12(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = ddc::get<DimR>(coord);
        const double p = ddc::get<DimP>(coord);
        return -r * std::sin(p);
    }

    double jacobian_21(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double p = ddc::get<DimP>(coord);
        return std::sin(p);
    }

    double jacobian_22(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = ddc::get<DimR>(coord);
        const double p = ddc::get<DimP>(coord);
        return r * std::cos(p);
    }


    void inv_jacobian_matrix(ddc::Coordinate<DimR, DimP> const& coord, Matrix_2x2& matrix)
            const final
    {
        const double r = ddc::get<DimR>(coord);
        const double p = ddc::get<DimP>(coord);
        assert(fabs(r) >= 1e-15);
        matrix[0][0] = std::cos(p);
        matrix[0][1] = std::sin(p);
        matrix[1][0] = -1 / r * std::sin(p);
        matrix[1][1] = 1 / r * std::cos(p);
    }

    double inv_jacobian_11(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double p = ddc::get<DimP>(coord);
        return std::cos(p);
    }

    double inv_jacobian_12(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double p = ddc::get<DimP>(coord);
        return std::sin(p);
    }

    double inv_jacobian_21(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = ddc::get<DimR>(coord);
        const double p = ddc::get<DimP>(coord);
        assert(fabs(r) >= 1e-15);
        return -1 / r * std::sin(p);
    }

    double inv_jacobian_22(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = ddc::get<DimR>(coord);
        const double p = ddc::get<DimP>(coord);
        assert(fabs(r) >= 1e-15);
        return 1 / r * std::cos(p);
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
    template <class Domain>
    void to_pseudo_cartesian_jacobian_center_matrix(Domain const& grid, Matrix_2x2& matrix) const
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
    template <class Domain>
    double to_pseudo_cartesian_jacobian_11_center(Domain const& grid) const
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
    template <class Domain>
    double to_pseudo_cartesian_jacobian_12_center(Domain const& grid) const
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
    template <class Domain>
    double to_pseudo_cartesian_jacobian_21_center(Domain const& grid) const
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
    template <class Domain>
    double to_pseudo_cartesian_jacobian_22_center(Domain const& grid) const
    {
        return 1.;
    }
};
