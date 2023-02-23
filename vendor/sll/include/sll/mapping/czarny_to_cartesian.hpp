#pragma once

#include <cmath>

#include <ddc/ddc.hpp>

#include <sll/mapping/curvilinear2d_to_cartesian.hpp>

template <class DimX, class DimY, class DimR, class DimP>
class CzarnyToCartesian : public Curvilinear2DToCartesian<DimX, DimY, DimR, DimP>
{
public:
    using cartesian_tag_x = DimX;
    using cartesian_tag_y = DimY;
    using circular_tag_r = DimR;
    using circular_tag_p = DimP;

private:
    double m_epsilon;
    double m_e;

public:
    CzarnyToCartesian(double epsilon, double e) : m_epsilon(epsilon), m_e(e) {}

    CzarnyToCartesian(CzarnyToCartesian const& other) = default;

    CzarnyToCartesian(CzarnyToCartesian&& x) = default;

    ~CzarnyToCartesian() = default;

    CzarnyToCartesian& operator=(CzarnyToCartesian const& x) = default;

    CzarnyToCartesian& operator=(CzarnyToCartesian&& x) = default;

    const double epsilon() const
    {
        return m_epsilon;
    }
    const double e() const
    {
        return m_e;
    }

    Coordinate<DimX, DimY> operator()(Coordinate<DimR, DimP> const& coord) const
    {
        const double r = get<DimR>(coord);
        const double theta = get<DimP>(coord);
        const double tmp1 = std::sqrt(m_epsilon * (m_epsilon + 2.0 * r * std::cos(theta)) + 1.0);

        const double x = (1.0 - tmp1) / m_epsilon;
        const double y = m_e * r * std::sin(theta)
                         / (std::sqrt(1.0 - 0.25 * m_epsilon * m_epsilon) * (2.0 - tmp1));

        return Coordinate<DimX, DimY>(x, y);
    }

    Coordinate<DimR, DimP> operator()(Coordinate<DimX, DimY> const& coord) const
    {
        const double x = get<DimX>(coord);
        const double y = get<DimY>(coord);
        const double ex = 1. + m_epsilon * x;
        const double ex2 = (m_epsilon * x * x - 2. * x - m_epsilon);
        const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
        const double xi = std::sqrt(xi2);
        const double r = std::sqrt(y * y * ex * ex / (m_e * m_e * xi2) + ex2 * ex2 * 0.25);
        const double theta
                = std::atan2(2. * y * ex, (m_e * xi * (m_epsilon * x * x - 2. * x - m_epsilon)));
        if (theta < 0)
            theta = 6.283185307179586 + theta;
        return Coordinate<DimR, DimP>(r, theta);
    }

    double jacobian_11(Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = get<DimR>(coord);
        const double theta = get<DimP>(coord);
        return -std::cos(theta)
               / std::sqrt(m_epsilon * (m_epsilon + 2.0 * r * std::cos(theta)) + 1.0);
    }

    double jacobian_12(Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = get<DimR>(coord);
        const double theta = get<DimP>(coord);
        return r * std::sin(theta)
               / std::sqrt(m_epsilon * (m_epsilon + 2.0 * r * std::cos(theta)) + 1.0);
    }

    double jacobian_21(Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = get<DimR>(coord);
        const double theta = get<DimP>(coord);

        const double sin_theta = std::sin(theta);
        const double cos_theta = std::cos(theta);
        const double xi2 = 1. / (1. - m_epsilon * m_epsilon * 0.25);
        const double xi = std::sqrt(xi2);
        const double tmp1 = std::sqrt(m_epsilon * (m_epsilon + 2.0 * r * cos_theta) + 1.0);
        const double tmp2 = 2.0 - tmp1;
        return m_e * m_epsilon * r * sin_theta * cos_theta * xi / (tmp2 * tmp2 * tmp1)
               + m_e * sin_theta * xi / tmp2;
    }

    double jacobian_22(Coordinate<DimR, DimP> const& coord) const final
    {
        const double r = get<DimR>(coord);
        const double theta = get<DimP>(coord);

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
