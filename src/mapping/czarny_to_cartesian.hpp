// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "tensor.hpp"
#include "view.hpp"

// Pre-declaration of analytical inverse
template <class X, class Y, class R, class Theta>
class CartesianToCzarny;

/**
 * @brief A class for describing the Czarny 2D mapping.
 *
 * The mapping @f$ (r,\theta)\mapsto (x,y) @f$ is defined by
 *
 * @f$ x(r,\theta) = \frac{1}{\epsilon} \left( 1 - \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))} \right),@f$
 *
 * @f$ y(r,\theta) = \frac{e\xi r \sin(\theta)}{2 -\sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))} },@f$
 *
 * with @f$ \xi = 1/\sqrt{1 - \epsilon^2 /4} @f$ and @f$ e @f$ and @f$ \epsilon @f$ given as parameters.
 * It and its Jacobian matrix are invertible everywhere except for @f$ r = 0 @f$.
 *
 * Its Jacobian coefficients are defined as follow
 *
 * @f$ J_{11}(r,\theta) = - \cos(\theta)\frac{1}{ \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))} } @f$
 *
 * @f$ J_{12}(r,\theta)  =  r\sin(\theta)\frac{1}{ \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))} } @f$
 *
 * @f$ J_{21}(r,\theta)  =  \cos(\theta)\frac{e\epsilon \xi r\sin(\theta)}{ \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))} \left(
 * 2 - \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))}  \right)^2 }
 * +  \sin(\theta)\frac{e\xi }{ 2- \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))} }@f$
 *
 * @f$ J_{22}(r,\theta)  =   r \sin(\theta)\frac{- e\epsilon \xi r \sin(\theta)}{ \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))} \left(
 * 2 - \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))} \right)^2 }
 * +  r\cos(\theta)\frac{e\xi }{ 2 -\sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))} }@f$.
 *
 *
 * and
 * @f$ \det(J(r, \theta)) = \frac{- r}{ \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))}}
 *  \frac{e\xi}{2 -  \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta))}}. @f$
 *
 */
template <class R, class Theta, class X, class Y>
class CzarnyToCartesian
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
    double m_epsilon;
    double m_e;

public:
    /**
     * @brief Instantiate a CzarnyToCartesian from parameters.
     *
     * @param[in] epsilon
     * 			The @f$ \epsilon @f$ parameter in the definition of the mapping CzarnyToCartesian.
     *
     * @param[in] e
     * 			The @f$ e @f$ parameter in the definition of the mapping CzarnyToCartesian.
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
    KOKKOS_FUNCTION CzarnyToCartesian(CzarnyToCartesian const& other)
        : m_epsilon(other.epsilon())
        , m_e(other.e())
    {
    }

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
    KOKKOS_FUNCTION double epsilon() const
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
    KOKKOS_FUNCTION double e() const
    {
        return m_e;
    }

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
        const double tmp1
                = Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * Kokkos::cos(theta)) + 1.0);

        const double x = (1.0 - tmp1) / m_epsilon;
        const double y = m_e * r * Kokkos::sin(theta)
                         / (Kokkos::sqrt(1.0 - 0.25 * m_epsilon * m_epsilon) * (2.0 - tmp1));

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
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        const double xi = Kokkos::sqrt(1. / (1. - m_epsilon * m_epsilon * 0.25));
        return -r / Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * Kokkos::cos(theta))) * m_e
               * xi
               / (2 - Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * Kokkos::cos(theta))));
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given independently with the function jacobian_component
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The Jacobian matrix.
     *
     */
    KOKKOS_FUNCTION DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix(
            Coord<R, Theta> const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);

        const double sin_theta = Kokkos::sin(theta);
        const double cos_theta = Kokkos::cos(theta);
        const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
        const double xi = Kokkos::sqrt(xi2);
        const double sqrt_eps = Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * cos_theta) + 1.0);
        const double sqrt_eps_2 = 2.0 - sqrt_eps;

        DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix;
        ddcHelper::get<X, R_cov>(jacobian_matrix) = -cos_theta / sqrt_eps;
        ddcHelper::get<X, Theta_cov>(jacobian_matrix) = r * sin_theta / sqrt_eps;
        ddcHelper::get<Y, R_cov>(jacobian_matrix) = m_e * m_epsilon * r * sin_theta * cos_theta * xi
                                                            / (sqrt_eps_2 * sqrt_eps_2 * sqrt_eps)
                                                    + m_e * sin_theta * xi / sqrt_eps_2;
        ddcHelper::get<Y, Theta_cov>(jacobian_matrix)
                = r
                  * (-m_e * m_epsilon * r * sin_theta * sin_theta * xi
                             / (sqrt_eps_2 * sqrt_eps_2 * sqrt_eps)
                     + m_e * cos_theta * xi / sqrt_eps_2);
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
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (1,1), i.e dx/dr
            return -Kokkos::cos(theta)
                   / Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * Kokkos::cos(theta)) + 1.0);
        } else if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, Theta_cov>) {
            // Component (1,2), i.e dx/dtheta
            return r * Kokkos::sin(theta)
                   / Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * Kokkos::cos(theta)) + 1.0);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (2,1), i.e dy/dr
            const double sin_theta = Kokkos::sin(theta);
            const double cos_theta = Kokkos::cos(theta);
            const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
            const double xi = Kokkos::sqrt(xi2);
            const double tmp1 = Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * cos_theta) + 1.0);
            const double tmp2 = 2.0 - tmp1;
            return m_e * m_epsilon * r * sin_theta * cos_theta * xi / (tmp2 * tmp2 * tmp1)
                   + m_e * sin_theta * xi / tmp2;
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, Theta_cov>) {
            // Component (2,2), i.e dy/dtheta
            const double sin_theta = Kokkos::sin(theta);
            const double cos_theta = Kokkos::cos(theta);
            const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
            const double xi = Kokkos::sqrt(xi2);
            const double tmp1 = Kokkos::sqrt(m_epsilon * (m_epsilon + 2.0 * r * cos_theta) + 1.0);
            const double tmp2 = 2.0 - tmp1;
            return r
                   * (-m_e * m_epsilon * r * sin_theta * sin_theta * xi / (tmp2 * tmp2 * tmp1)
                      + m_e * cos_theta * xi / tmp2);
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

        assert(r >= 1e-15);

        const double sin_theta = Kokkos::sin(theta);
        const double cos_theta = Kokkos::cos(theta);
        const double xi = Kokkos::sqrt(1. / (1. - m_epsilon * m_epsilon * 0.25));
        const double divisor = 2 - Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * cos_theta));

        const double fact_1 = 1 / Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * cos_theta));
        const double fact_2 = m_e * m_epsilon * xi * r * sin_theta * fact_1 / divisor / divisor;
        const double fact_3 = m_e * xi / divisor;

        DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X_cov, Y_cov>> matrix;
        ddcHelper::get<R, X_cov>(matrix)
                = -1 / fact_1 * (-sin_theta * fact_2 + cos_theta * fact_3) / fact_3;
        ddcHelper::get<R, Y_cov>(matrix) = sin_theta / fact_3;
        ddcHelper::get<Theta, X_cov>(matrix)
                = 1 / r / fact_1 * (cos_theta * fact_2 + sin_theta * fact_3) / fact_3;
        ddcHelper::get<Theta, Y_cov>(matrix) = 1 / r * cos_theta / fact_3;
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

        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        const double cos_theta = Kokkos::cos(theta);
        const double xi = Kokkos::sqrt(1. / (1. - m_epsilon * m_epsilon * 0.25));
        const double divisor = 2 - Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * cos_theta));

        if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, X_cov>) {
            //Compute the (1,1) coefficient of the inverse Jacobian matrix
            const double sin_theta = Kokkos::sin(theta);

            const double fact_1
                    = 1 / Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * cos_theta));
            const double fact_2 = m_e * m_epsilon * xi * r * sin_theta * fact_1 / divisor / divisor;
            const double fact_3 = m_e * xi / divisor;

            return -1 / fact_1 * (-sin_theta * fact_2 + cos_theta * fact_3) / fact_3;
        } else if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, Y_cov>) {
            //Compute the (1,2) coefficient of the Jacobian matrix, i.e J^x_theta.
            const double sin_theta = Kokkos::sin(theta);

            const double fact_3 = m_e * xi / divisor;
            return sin_theta / fact_3;
        } else if constexpr (std::is_same_v<IndexTag1, Theta> && std::is_same_v<IndexTag2, X_cov>) {
            //Compute the (2,1) coefficient of the inverse Jacobian matrix.
            assert(r >= 1e-15);

            const double sin_theta = Kokkos::sin(theta);

            const double fact_1
                    = 1 / Kokkos::sqrt(1 + m_epsilon * (m_epsilon + 2.0 * r * cos_theta));
            const double fact_2 = m_e * m_epsilon * xi * r * sin_theta * fact_1 / divisor / divisor;
            const double fact_3 = m_e * xi / divisor;

            return 1 / r / fact_1 * (cos_theta * fact_2 + sin_theta * fact_3) / fact_3;
        } else {
            //Compute the (2,2) coefficient of the inverse Jacobian matrix.
            assert(r >= 1e-15);

            const double cos_theta = Kokkos::cos(theta);

            const double fact_3 = m_e * xi / divisor;
            return 1 / r * cos_theta / fact_3;
        }
    }


    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
    CartesianToCzarny<X, Y, R, Theta> get_inverse_mapping() const
    {
        return CartesianToCzarny<X, Y, R, Theta>(m_epsilon, m_e);
    }
};

namespace mapping_detail {
template <class X, class Y, class R, class Theta, class ExecSpace>
struct MappingAccessibility<ExecSpace, CzarnyToCartesian<R, Theta, X, Y>> : std::true_type
{
};

template <class X, class Y, class R, class Theta>
struct IsCurvilinear2DMapping<CzarnyToCartesian<X, Y, R, Theta>> : std::true_type
{
};

template <class X, class Y, class R, class Theta>
struct SingularOPointInvJacobian<CzarnyToCartesian<R, Theta, X, Y>> : std::true_type
{
};

} // namespace mapping_detail
