#pragma once

#include <cmath>

#include <ddc/ddc.hpp>

#include "analytical_invertible_curvilinear2d_to_cartesian.hpp"


/**
 * @brief A class for describing the Czarny 2D mapping.
 *
 * The mapping @f$ (r,\theta)\mapsto (x,y) @f$ is defined by
 *
 * @f$ x(r,\theta) = \frac{1}{\epsilon} \left( 1 - \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} \right),@f$
 *
 * @f$ y(r,\theta) = \frac{e\xi r \sin(\theta)}{2 -\sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} },@f$
 *
 * with @f$ \xi = 1/\sqrt{1 - \epsilon^2 /4} @f$ and @f$ e @f$ and @f$ \epsilon @f$ given as parameters.
 * It and its Jacobian matrix are invertible everywhere except for @f$ r = 0 @f$.
 *
 * Its Jacobian coefficients are defined as follow
 *
 * @f$ J_{11}(r,\theta) = - \frac{\cos(\theta)}{ \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} } @f$
 *
 * @f$ J_{12}(r,\theta)  =  \frac{r\sin(\theta)}{ \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} } @f$
 *
 * @f$ J_{21}(r,\theta)  =  \frac{e\epsilon \xi r\sin(\theta)\cos(\theta)}{ \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} \left(
 * 2 - \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)}  \right)^2 }
 * +  \frac{e\xi \sin(\theta)}{ \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} }@f$
 *
 * @f$ J_{22}(r,\theta)  =   \frac{e\epsilon \xi r^2 \sin^2(\theta)}{ \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} \left(
 * 2 - \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} \right)^2 }
 * +  \frac{e\xi r\cos(\theta)}{ \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} }@f$.
 *
 *
 * @see AnalyticalInvertibleCurvilinear2DToCartesian
 */
template <class DimX, class DimY, class DimR, class DimP>
class CzarnyToCartesian
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

private:
    double m_epsilon;
    double m_e;

public:
    /**
     * @brief Instantiate a CzarnyToCartesian from parameters.
     *
     * @param[in] epsilon
     * 			The @f$ \epsilon @f$ parameter in the definition of the mapping 'CzarnyToCartesian'.
     *
     * @param[in] e
     * 			The @f$ e @f$ parameter in the definition of the mapping 'CzarnyToCartesian'.
     *
     * @see CzarnyToCartesian
     */
    CzarnyToCartesian(double epsilon, double e) : m_epsilon(epsilon), m_e(e) {}

    /**
     * @brief Instantiate a CzarnyToCartesian from another CzarnyToCartesian (lvalue).
     *
     * @param[in] other
     * 		CzarnyToCartesian mapping used to instantiate the new one.
     */
    CzarnyToCartesian(CzarnyToCartesian const& other) = default;

    /**
     * @brief Instantiate a CzarnyToCartesian from another temporary CzarnyToCartesian (rvalue).
     *
     * @param[in] x
     * 		CzarnyToCartesian mapping used to instantiate the new one.
     */
    CzarnyToCartesian(CzarnyToCartesian&& x) = default;

    ~CzarnyToCartesian() = default;

    /**
     * @brief Assign a CzarnyToCartesian from another CzarnyToCartesian (lvalue).
     *
     * @param[in] x
     * 		CzarnyToCartesian mapping used to assign.
     *
     * @return The CzarnyToCartesian assigned.
     */
    CzarnyToCartesian& operator=(CzarnyToCartesian const& x) = default;

    /**
     * @brief Assign a CzarnyToCartesian from another temporary CzarnyToCartesian (rvalue).
     *
     * @param[in] x
     * 		CzarnyToCartesian mapping used to assign.
     *
     * @return The CzarnyToCartesian assigned.
     */
    CzarnyToCartesian& operator=(CzarnyToCartesian&& x) = default;

    /**
     * @brief Return the @f$ \epsilon @f$ parameter.
     *
     * @return The value of @f$ \epsilon @f$.
     *
     * @see CzarnyToCartesian
     */
    const double epsilon() const
    {
        return m_epsilon;
    }

    /**
     * @brief Return the @f$ e @f$ parameter.
     *
     * @return The value of @f$ e @f$.
     *
     * @see CzarnyToCartesian
     */
    const double e() const
    {
        return m_e;
    }

    ddc::Coordinate<DimX, DimY> operator()(ddc::Coordinate<DimR, DimP> const& coord) const
    {
        const double r = ddc::get<DimR>(coord);
        const double theta = ddc::get<DimP>(coord);
        const double tmp1 = std::sqrt(m_epsilon * (m_epsilon + 2.0 * r * std::cos(theta)) + 1.0);

        const double x = (1.0 - tmp1) / m_epsilon;
        const double y = m_e * r * std::sin(theta)
                         / (std::sqrt(1.0 - 0.25 * m_epsilon * m_epsilon) * (2.0 - tmp1));

        return ddc::Coordinate<DimX, DimY>(x, y);
    }

    ddc::Coordinate<DimR, DimP> operator()(ddc::Coordinate<DimX, DimY> const& coord) const
    {
        const double x = ddc::get<DimX>(coord);
        const double y = ddc::get<DimY>(coord);
        const double ex = 1. + m_epsilon * x;
        const double ex2 = (m_epsilon * x * x - 2. * x - m_epsilon);
        const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
        const double xi = std::sqrt(xi2);
        const double r = std::sqrt(y * y * ex * ex / (m_e * m_e * xi2) + ex2 * ex2 * 0.25);
        double theta
                = std::atan2(2. * y * ex, (m_e * xi * (m_epsilon * x * x - 2. * x - m_epsilon)));
        if (theta < 0) {
            theta = 2 * M_PI + theta;
        }
        return ddc::Coordinate<DimR, DimP>(r, theta);
    }

    void jacobian_matrix(ddc::Coordinate<DimR, DimP> const& coord, Matrix_2x2& matrix) const final
    {
        const double r = ddc::get<DimR>(coord);
        const double theta = ddc::get<DimP>(coord);

        const double sin_theta = std::sin(theta);
        const double cos_theta = std::cos(theta);
        const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
        const double xi = std::sqrt(xi2);
        const double sqrt_eps = std::sqrt(m_epsilon * (m_epsilon + 2.0 * r * cos_theta) + 1.0);
        const double sqrt_eps_2 = 2.0 - sqrt_eps;

        matrix[0][0] = -cos_theta / sqrt_eps;
        matrix[0][1] = r * sin_theta / sqrt_eps;
        matrix[1][0] = m_e * m_epsilon * r * sin_theta * cos_theta * xi
                               / (sqrt_eps_2 * sqrt_eps_2 * sqrt_eps)
                       + m_e * sin_theta * xi / sqrt_eps_2;
        matrix[1][1] = r
                       * (-m_e * m_epsilon * r * sin_theta * sin_theta * xi
                                  / (sqrt_eps_2 * sqrt_eps_2 * sqrt_eps)
                          + m_e * cos_theta * xi / sqrt_eps_2);
    }

    double jacobian_11(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = ddc::get<DimR>(coord);
        const double theta = ddc::get<DimP>(coord);
        return -std::cos(theta)
               / std::sqrt(m_epsilon * (m_epsilon + 2.0 * r * std::cos(theta)) + 1.0);
    }

    double jacobian_12(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = ddc::get<DimR>(coord);
        const double theta = ddc::get<DimP>(coord);
        return r * std::sin(theta)
               / std::sqrt(m_epsilon * (m_epsilon + 2.0 * r * std::cos(theta)) + 1.0);
    }

    double jacobian_21(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = ddc::get<DimR>(coord);
        const double theta = ddc::get<DimP>(coord);

        const double sin_theta = std::sin(theta);
        const double cos_theta = std::cos(theta);
        const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
        const double xi = std::sqrt(xi2);
        const double tmp1 = std::sqrt(m_epsilon * (m_epsilon + 2.0 * r * cos_theta) + 1.0);
        const double tmp2 = 2.0 - tmp1;
        return m_e * m_epsilon * r * sin_theta * cos_theta * xi / (tmp2 * tmp2 * tmp1)
               + m_e * sin_theta * xi / tmp2;
    }

    double jacobian_22(ddc::Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = ddc::get<DimR>(coord);
        const double theta = ddc::get<DimP>(coord);

        const double sin_theta = std::sin(theta);
        const double cos_theta = std::cos(theta);
        const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
        const double xi = std::sqrt(xi2);
        const double tmp1 = std::sqrt(m_epsilon * (m_epsilon + 2.0 * r * cos_theta) + 1.0);
        const double tmp2 = 2.0 - tmp1;
        return r
               * (-m_e * m_epsilon * r * sin_theta * sin_theta * xi / (tmp2 * tmp2 * tmp1)
                  + m_e * cos_theta * xi / tmp2);
    }
};
