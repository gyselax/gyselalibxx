// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <cassert>

#include <sll/view.hpp>

/**
 * An operator to calculate the Jacobian matrix and its inverse.
 * All operators which can calculate terms of the Jacobian matrix should inherit from this class.
 *
 * @tparam PositionCoordinate The type of the coordinate at which the Jacobian matrix can be calculated.
 */
template <class PositionCoordinate>
class Jacobian
{
public:
    /// The type of the Jacobian matrix and its inverse
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

public:
    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     * 			The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    virtual double jacobian(PositionCoordinate const& coord) const = 0;

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
     *
     *
     * @see Jacobian::jacobian_11
     * @see Jacobian::jacobian_12
     * @see Jacobian::jacobian_21
     * @see Jacobian::jacobian_22
     */
    virtual void jacobian_matrix(PositionCoordinate const& coord, Matrix_2x2& matrix) const = 0;

    /**
     * @brief Compute the (1,1) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (1,1) coefficient of the Jacobian matrix is given by @f$ \frac{\partial x}{\partial r} @f$.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (1,1) coefficient of the Jacobian matrix.
     */
    virtual double jacobian_11(PositionCoordinate const& coord) const = 0;

    /**
     * @brief Compute the (1,2) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (1,2) coefficient of the Jacobian matrix is given by @f$ \frac{\partial x}{\partial \theta} @f$.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (1,2) coefficient of the Jacobian matrix.
     */
    virtual double jacobian_12(PositionCoordinate const& coord) const = 0;

    /**
     * @brief Compute the (2,1) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (2,1) coefficient of the Jacobian matrix is given by @f$ \frac{\partial y}{\partial r} @f$.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix. .
     *
     * @return A double with the value of the (2,1) coefficient of the Jacobian matrix.
     */
    virtual double jacobian_21(PositionCoordinate const& coord) const = 0;

    /**
     * @brief Compute the (2,2) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (2,2) coefficient of the Jacobian matrix is given by @f$ \frac{\partial y}{\partial \theta} @f$.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (2,2) coefficient of the Jacobian matrix.
     */
    virtual double jacobian_22(PositionCoordinate const& coord) const = 0;

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
    virtual void inv_jacobian_matrix(PositionCoordinate const& coord, Matrix_2x2& matrix) const = 0;

    /**
     * @brief Compute the (1,1) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (1,1) coefficient of the inverse Jacobian matrix.
     */
    virtual double inv_jacobian_11(PositionCoordinate const& coord) const = 0;

    /**
     * @brief Compute the (1,2) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (1,2) coefficient of the inverse Jacobian matrix.
     */
    virtual double inv_jacobian_12(PositionCoordinate const& coord) const = 0;

    /**
     * @brief Compute the (2,1) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (2,1) coefficient of the inverse Jacobian matrix.
     */
    virtual double inv_jacobian_21(PositionCoordinate const& coord) const = 0;

    /**
     * @brief Compute the (2,2) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (2,2) coefficient of the inverse Jacobian matrix.
     */
    virtual double inv_jacobian_22(PositionCoordinate const& coord) const = 0;
};

/**
 * A specialisation of Jacobian to handle non-analytical terms. In this case the inverse and the determinant
 * are calculated from the Jacobian matrix in the same way regardless of the implementation of the calculation
 * of the Jacobian itself.
 *
 * @tparam PositionCoordinate The type of the coordinate at which the Jacobian matrix can be calculated.
 */
template <class PositionCoordinate>
class NonAnalyticalJacobian : public Jacobian<PositionCoordinate>
{
public:
    /// The type of the Jacobian matrix and its inverse
    using Matrix_2x2 = typename Jacobian<PositionCoordinate>::Matrix_2x2;

public:
    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     * 			The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    double jacobian(PositionCoordinate const& coord) const final
    {
        const double j_rr = this->jacobian_11(coord);
        const double j_rp = this->jacobian_12(coord);
        const double j_pr = this->jacobian_21(coord);
        const double j_pp = this->jacobian_22(coord);
        return j_rr * j_pp - j_rp * j_pr;
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
    void inv_jacobian_matrix(PositionCoordinate const& coord, Matrix_2x2& matrix) const final
    {
        double jacob = jacobian(coord);
        assert(fabs(jacob) > 1e-15);
        matrix[0][0] = this->jacobian_22(coord) / jacob;
        matrix[0][1] = -this->jacobian_12(coord) / jacob;
        matrix[1][0] = -this->jacobian_21(coord) / jacob;
        matrix[1][1] = this->jacobian_11(coord) / jacob;
    }

    /**
     * @brief Compute the (1,1) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (1,1) coefficient of the inverse Jacobian matrix.
     */
    double inv_jacobian_11(PositionCoordinate const& coord) const final
    {
        double jacob = jacobian(coord);
        assert(fabs(jacob) > 1e-15);
        return this->jacobian_22(coord) / jacob;
    }

    /**
     * @brief Compute the (1,2) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (1,2) coefficient of the inverse Jacobian matrix.
     */
    double inv_jacobian_12(PositionCoordinate const& coord) const final
    {
        double jacob = jacobian(coord);
        assert(fabs(jacob) > 1e-15);
        return -this->jacobian_12(coord) / jacob;
    }

    /**
     * @brief Compute the (2,1) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (2,1) coefficient of the inverse Jacobian matrix.
     */
    double inv_jacobian_21(PositionCoordinate const& coord) const final
    {
        double jacob = jacobian(coord);
        assert(fabs(jacob) > 1e-15);
        return -this->jacobian_21(coord) / jacob;
    }

    /**
     * @brief Compute the (2,2) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (2,2) coefficient of the inverse Jacobian matrix.
     */
    double inv_jacobian_22(PositionCoordinate const& coord) const final
    {
        double jacob = jacobian(coord);
        assert(fabs(jacob) > 1e-15);
        return this->jacobian_11(coord) / jacob;
    }
};
