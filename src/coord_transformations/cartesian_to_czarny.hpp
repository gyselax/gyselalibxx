// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_aliases.hpp"
#include "view.hpp"

// Pre-declaration of analytical inverse
template <class R, class Theta, class X, class Y>
class CzarnyToCartesian;

/**
 * @brief A class for describing the Czarny 2D mapping.
 *
 * The mapping @f$ (x,y)  \mapsto (r,\theta)@f$ is defined by
 *
 * @f$ r(x,y) = \sqrt{\frac{\hat{y}^2 (1+\epsilon \hat{x})^2}{e^2\xi^2+0.25(\epsilon \hat{x}^2-2\hat{x}-\epsilon)^2}},@f$
 *
 * @f$ \theta (x,y)) = atan2(2. \hat{y} (1+\epsilon \hat{x}), (e \xi (\epsilon \hat{x}^2 - 2\hat{x}-\epsilon))), @f$
 *
 * with @f$ \hat{x} = x-x_0 @f$, @f$ \hat{y} = y-y_0 @f$, @f$ \xi = 1/\sqrt{1 - \epsilon^2 /4} @f$ and @f$ e @f$ and @f$ \epsilon @f$ given as parameters.
 *
 * See O. Czarny and G. Huysmans. Bézier surfaces and finite elements for MHD simulations.
 * Journal of Computational Physics, 227(16):7423–7445, 2008. doi:10.1016/j.jcp.2008.04.001.
 */
template <class X, class Y, class R, class Theta>
class CartesianToCzarny
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

private:
    double m_epsilon;
    double m_e;
    Coord<X, Y> m_o_point;

public:
    /**
     * @brief Instantiate a CartesianToCzarny from parameters.
     *
     * @param[in] epsilon
     * 			The @f$ \epsilon @f$ parameter in the definition of the mapping CartesianToCzarny.
     * @param[in] e
     * 			The @f$ e @f$ parameter in the definition of the mapping CartesianToCzarny.
     * @param[in] o_point The (x,y)-coordinate of the O-point ((0,0) by default).
     *
     * @see CartesianToCzarny
     */
    explicit KOKKOS_FUNCTION CartesianToCzarny(
            double epsilon,
            double e,
            Coord<X, Y> o_point = Coord<X, Y>(0.0, 0.0))
        : m_epsilon(epsilon)
        , m_e(e)
        , m_o_point(o_point)
    {
    }

    /**
     * @brief Instantiate a CartesianToCzarny from another CartesianToCzarny (lvalue).
     *
     * @param[in] other
     * 		CartesianToCzarny mapping used to instantiate the new one.
     */
    KOKKOS_DEFAULTED_FUNCTION CartesianToCzarny(CartesianToCzarny const& other) = default;

    /**
     * @brief Instantiate a CartesianToCzarny from another temporary CartesianToCzarny (rvalue).
     *
     * @param[in] x
     * 		CartesianToCzarny mapping used to instantiate the new one.
     */
    CartesianToCzarny(CartesianToCzarny&& x) = default;

    KOKKOS_DEFAULTED_FUNCTION ~CartesianToCzarny() = default;

    /**
     * @brief Assign a CartesianToCzarny from another CartesianToCzarny (lvalue).
     *
     * @param[in] x
     * 		CartesianToCzarny mapping used to assign.
     *
     * @return The CartesianToCzarny assigned.
     */
    CartesianToCzarny& operator=(CartesianToCzarny const& x) = default;

    /**
     * @brief Assign a CartesianToCzarny from another temporary CartesianToCzarny (rvalue).
     *
     * @param[in] x
     * 		CartesianToCzarny mapping used to assign.
     *
     * @return The CartesianToCzarny assigned.
     */
    CartesianToCzarny& operator=(CartesianToCzarny&& x) = default;

    /**
     * @brief Return the @f$ \epsilon @f$ parameter.
     *
     * @return The value of @f$ \epsilon @f$.
     *
     * @see CartesianToCzarny
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
     * @see CartesianToCzarny
     */
    KOKKOS_FUNCTION double e() const
    {
        return m_e;
    }

    /**
     * @brief Convert the coordinate (x,y) to the equivalent @f$ (r, \theta) @f$ coordinate.
     *
     * @param[in] coord The coordinate to be converted.
     *
     * @return The equivalent coordinate.
     */
    KOKKOS_FUNCTION Coord<R, Theta> operator()(Coord<X, Y> const& coord) const
    {
        const double x = ddc::get<X>(coord) - ddc::get<X>(m_o_point);
        const double y = ddc::get<Y>(coord) - ddc::get<Y>(m_o_point);
        const double ex = 1. + m_epsilon * x;
        const double ex2 = (m_epsilon * x * x - 2. * x - m_epsilon);
        const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
        const double xi = Kokkos::sqrt(xi2);
        const double r = Kokkos::sqrt(y * y * ex * ex / (m_e * m_e * xi2) + ex2 * ex2 * 0.25);
        double theta = Kokkos::atan2(2. * y * ex, (m_e * xi * ex2));
        if (theta < 0) {
            theta = 2 * M_PI + theta;
        }
        return Coord<R, Theta>(r, theta);
    }

    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
    KOKKOS_INLINE_FUNCTION CzarnyToCartesian<R, Theta, X, Y> get_inverse_mapping() const
    {
        return CzarnyToCartesian<R, Theta, X, Y>(m_epsilon, m_e, m_o_point);
    }
};

namespace mapping_detail {
template <class X, class Y, class R, class Theta, class ExecSpace>
struct MappingAccessibility<ExecSpace, CartesianToCzarny<X, Y, R, Theta>> : std::true_type
{
};
} // namespace mapping_detail
