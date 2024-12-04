// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include <sll/view.hpp>

#include "mapping_tools.hpp"


/**
 * @brief A class for describing the Czarny 2D mapping.
 *
 * The mapping @f$ (x,y)  \mapsto (r,\theta)@f$ is defined by
 *
 * @f$ r(x,y) = \sqrt{\frac{y^2 (1+\epsilon x)^2}{e^2\xi^2+0.25(\epsilon x^2-2x-\epsilon)^2}},@f$
 *
 * @f$ \theta (x,y)) = atan2(2. y (1+\epsilon x), (e \xi (\epsilon x^2 - 2x-\epsilon))), @f$
 *
 * with @f$ \xi = 1/\sqrt{1 - \epsilon^2 /4} @f$ and @f$ e @f$ and @f$ \epsilon @f$ given as parameters.
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
    using CoordArg = ddc::Coordinate<X, Y>;
    /// The type of the result of the function described by this mapping
    using CoordResult = ddc::Coordinate<R, Theta>;

private:
    double m_epsilon;
    double m_e;

public:
    /**
     * @brief Instantiate a CartesianToCzarny from parameters.
     *
     * @param[in] epsilon
     * 			The @f$ \epsilon @f$ parameter in the definition of the mapping CartesianToCzarny.
     *
     * @param[in] e
     * 			The @f$ e @f$ parameter in the definition of the mapping CartesianToCzarny.
     *
     * @see CartesianToCzarny
     */
    CartesianToCzarny(double epsilon, double e) : m_epsilon(epsilon), m_e(e) {}

    /**
     * @brief Instantiate a CartesianToCzarny from another CartesianToCzarny (lvalue).
     *
     * @param[in] other
     * 		CartesianToCzarny mapping used to instantiate the new one.
     */
    KOKKOS_FUNCTION CartesianToCzarny(CartesianToCzarny const& other)
        : m_epsilon(other.epsilon())
        , m_e(other.e())
    {
    }

    /**
     * @brief Instantiate a CartesianToCzarny from another temporary CartesianToCzarny (rvalue).
     *
     * @param[in] x
     * 		CartesianToCzarny mapping used to instantiate the new one.
     */
    CartesianToCzarny(CartesianToCzarny&& x) = default;

    ~CartesianToCzarny() = default;

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
    KOKKOS_FUNCTION ddc::Coordinate<R, Theta> operator()(ddc::Coordinate<X, Y> const& coord) const
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        const double ex = 1. + m_epsilon * x;
        const double ex2 = (m_epsilon * x * x - 2. * x - m_epsilon);
        const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
        const double xi = Kokkos::sqrt(xi2);
        const double r = Kokkos::sqrt(y * y * ex * ex / (m_e * m_e * xi2) + ex2 * ex2 * 0.25);
        double theta = Kokkos::atan2(2. * y * ex, (m_e * xi * ex2));
        if (theta < 0) {
            theta = 2 * M_PI + theta;
        }
        return ddc::Coordinate<R, Theta>(r, theta);
    }
};

namespace mapping_detail {
template <class X, class Y, class R, class Theta, class ExecSpace>
struct MappingAccessibility<ExecSpace, CartesianToCzarny<X, Y, R, Theta>> : std::true_type
{
};
} // namespace mapping_detail
