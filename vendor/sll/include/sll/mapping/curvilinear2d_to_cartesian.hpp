#pragma once

#include <array>
#include <cassert>

#include <ddc/ddc.hpp>



/**
 * @brief A class for describing curvilinear 2D mappings from the logical domain to the physical domain.
 * */
template <class DimX, class DimY, class DimR, class DimP>
class Curvilinear2DToCartesian
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
    using curvilinear_tag_r = DimR;
    /**
     * @brief Indicate the second logical coordinate.
     */
    using curvilinear_tag_p = DimP;
    /**
     * @brief Define a 2x2 matrix with an 2D array of an 2D array.
     */
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

public:
    Curvilinear2DToCartesian() = default;

    /**
     * @brief Instantiate a Curvilinear2DToCartesian from another
     * Curvilinear2DToCartesian (lvalue).
     *
     * @param[in] other
     * 		Curvilinear2DToCartesian mapping used to instantiate the new one.
     */
    Curvilinear2DToCartesian(Curvilinear2DToCartesian const& other) = default;

    /**
     * @brief Instantiate a Curvilinear2DToCartesian from another temporary
     *  Curvilinear2DToCartesian (rvalue).
     *
     * The input x
     *
     * @param[in] x
     * 		Curvilinear2DToCartesian mapping used to instantiate the new one.
     */
    Curvilinear2DToCartesian(Curvilinear2DToCartesian&& x) = default;

    virtual ~Curvilinear2DToCartesian() = default;

    /**
     * @brief Assign a Curvilinear2DToCartesian from another Curvilinear2DToCartesian (lvalue).
     *
     * @param[in] x
     * 		Curvilinear2DToCartesian mapping used to assign.
     *
     * @return The Curvilinear2DToCartesian assigned.
     */
    Curvilinear2DToCartesian& operator=(Curvilinear2DToCartesian const& x) = default;

    /**
     * @brief Assign a Curvilinear2DToCartesian from another temporary Curvilinear2DToCartesian (rvalue).
     *
     * @param[in] x
     * 		Curvilinear2DToCartesian mapping used to assign.
     *
     * @return The Curvilinear2DToCartesian assigned.
     */
    Curvilinear2DToCartesian& operator=(Curvilinear2DToCartesian&& x) = default;

    /**
     * @brief Compute the physical coordinates from the logical coordinates.
     *
     * The mapping from the logical domain @f$ (r,\theta) @f$ to the physical domain @f$ (x,y) @f$
     * must always be well defined.
     * However the inverse mapping is not always well defined especially at the center point
     * or can be costly to inverse.
     *
     * @param[in] coord
     * 			The coordinates in the logical domain.
     *
     * @return The coordinates in the physical domain.
     *
     */
    virtual ddc::Coordinate<DimX, DimY> operator()(
            ddc::Coordinate<DimR, DimP> const& coord) const = 0;

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     * 			The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    virtual double jacobian(ddc::Coordinate<DimR, DimP> const& coord) const
    {
        const double j_rr = jacobian_11(coord);
        const double j_rp = jacobian_12(coord);
        const double j_pr = jacobian_21(coord);
        const double j_pp = jacobian_22(coord);
        return j_rr * j_pp - j_rp * j_pr;
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given indendently with the functions
     * 'jacobian_11', 'jacobian_12', ' jacobian_21' and 'jacobian_22'.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @param[out] matrix
     * 				The Jacobian matrix returned.
     *
     *
     * @see Curvilinear2DToCartesian::jacobian_11
     * @see Curvilinear2DToCartesian::jacobian_12
     * @see Curvilinear2DToCartesian::jacobian_21
     * @see Curvilinear2DToCartesian::jacobian_22
     */
    virtual void jacobian_matrix(ddc::Coordinate<DimR, DimP> const& coord, Matrix_2x2& matrix)
            const = 0;

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
    virtual double jacobian_11(ddc::Coordinate<DimR, DimP> const& coord) const = 0;
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
    virtual double jacobian_12(ddc::Coordinate<DimR, DimP> const& coord) const = 0;
    /**
     * @brief Compute the (2,1) coefficient of the Jacobian matrix.
     *
     *For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (2,1) coefficient of the Jacobian matrix is given by @f$ \frac{\partial y}{\partial r} @f$.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix. .
     *
     * @return A double with the value of the (2,1) coefficient of the Jacobian matrix.
     */
    virtual double jacobian_21(ddc::Coordinate<DimR, DimP> const& coord) const = 0;
    /**
     * @brief Compute the (2,2) coefficient of the Jacobian matrix.
     *
     *For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (2,2) coefficient of the Jacobian matrix is given by @f$ \frac{\partial y}{\partial \theta} @f$.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (2,2) coefficient of the Jacobian matrix.
     */
    virtual double jacobian_22(ddc::Coordinate<DimR, DimP> const& coord) const = 0;

    /**
     * @brief Compute full inverse Jacobian matrix.
     *
     * For some computations, we need the complete inverse Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given indendently with the functions
     * 'inv_jacobian_11', 'inv_jacobian_12', 'inv_jacobian_21' and 'inv_jacobian_22'.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @param[out] matrix
     * 				The inverse Jacobian matrix returned.
     *
     *
     * @see Curvilinear2DToCartesian::inv_jacobian_11
     * @see Curvilinear2DToCartesian::inv_jacobian_12
     * @see Curvilinear2DToCartesian::inv_jacobian_21
     * @see Curvilinear2DToCartesian::inv_jacobian_22
     */
    virtual void inv_jacobian_matrix(ddc::Coordinate<DimR, DimP> const& coord, Matrix_2x2& matrix)
            const
    {
        double jacob = jacobian(coord);
        assert(fabs(jacobian(coord)) >= 1e-15);
        matrix[0][0] = jacobian_22(coord) / jacob;
        matrix[0][1] = -jacobian_12(coord) / jacob;
        matrix[1][0] = -jacobian_21(coord) / jacob;
        matrix[1][1] = jacobian_11(coord) / jacob;
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
    virtual double inv_jacobian_11(ddc::Coordinate<DimR, DimP> const& coord) const
    {
        assert(fabs(jacobian(coord)) >= 1e-15);
        return jacobian_22(coord) / jacobian(coord);
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
    virtual double inv_jacobian_12(ddc::Coordinate<DimR, DimP> const& coord) const
    {
        assert(fabs(jacobian(coord)) >= 1e-15);
        return -jacobian_12(coord) / jacobian(coord);
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
    virtual double inv_jacobian_21(ddc::Coordinate<DimR, DimP> const& coord) const
    {
        assert(fabs(jacobian(coord)) >= 1e-15);
        return -jacobian_21(coord) / jacobian(coord);
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
    virtual double inv_jacobian_22(ddc::Coordinate<DimR, DimP> const& coord) const
    {
        assert(fabs(jacobian(coord)) >= 1e-15);
        return jacobian_11(coord) / jacobian(coord);
    }


    /**
     * @brief Compute the metric tensor assignd to the mapping.
     *
     * The metric tensor matrix is defined for mapping whose the Jacobian matrix is called
     * @f$ J_{\mathcal{F}} @f$ as @f$ G = (J_{\mathcal{F}})^T J_{\mathcal{F}} @f$.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the metric tensor.
     * @param[out] matrix
     * 				The metric tensor matrix.
     */
    virtual void metric_tensor(ddc::Coordinate<DimR, DimP> const& coord, Matrix_2x2& matrix) const
    {
        const double J_rr = jacobian_11(coord);
        const double J_rp = jacobian_12(coord);
        const double J_pr = jacobian_21(coord);
        const double J_pp = jacobian_22(coord);
        matrix[0][0] = (J_rr * J_rr + J_pr * J_pr);
        matrix[0][1] = (J_rr * J_rp + J_pr * J_pp);
        matrix[1][0] = (J_rr * J_rp + J_pr * J_pp);
        matrix[1][1] = (J_rp * J_rp + J_pp * J_pp);
    }

    /**
     * @brief Compute the inverse metric tensor associated to the mapping.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the metric tensor.
     * @param[out] matrix
     * 				The metric tensor matrix.
     */
    virtual void inverse_metric_tensor(ddc::Coordinate<DimR, DimP> const& coord, Matrix_2x2& matrix)
            const
    {
        assert(fabs(ddc::get<DimR>(coord)) >= 1e-15);
        const double J_rr = jacobian_11(coord);
        const double J_rp = jacobian_12(coord);
        const double J_pr = jacobian_21(coord);
        const double J_pp = jacobian_22(coord);
        const double jacob_2 = jacobian(coord) * jacobian(coord);
        matrix[0][0] = (J_rp * J_rp + J_pp * J_pp) / jacob_2;
        matrix[0][1] = (-J_rr * J_rp - J_pr * J_pp) / jacob_2;
        matrix[1][0] = (-J_rr * J_rp - J_pr * J_pp) / jacob_2;
        matrix[1][1] = (J_rr * J_rr + J_pr * J_pr) / jacob_2;
    }

    /**
     * @brief Compute the covariant vector from the contravariant vector
     *
     * @param[in] contravariant_vector
     * 				The metric tensor matrix.
     * @param[in] coord
     * 				The coordinate where we want to compute the convariant vector.
     *
     * @return A vector of the covariant
     */
    std::array<double, 2> to_covariant(
            std::array<double, 2> const& contravariant_vector,
            ddc::Coordinate<DimR, DimP> const& coord) const
    {
        Matrix_2x2 inv_metric_tensor;
        inverse_metric_tensor(coord, inv_metric_tensor);
        std::array<double, 2> covariant_vector;
        covariant_vector[0] = inv_metric_tensor[0][0] * contravariant_vector[0]
                              + inv_metric_tensor[0][1] * contravariant_vector[1];
        covariant_vector[1] = inv_metric_tensor[1][0] * contravariant_vector[0]
                              + inv_metric_tensor[1][1] * contravariant_vector[1];
        return covariant_vector;
    }
};
