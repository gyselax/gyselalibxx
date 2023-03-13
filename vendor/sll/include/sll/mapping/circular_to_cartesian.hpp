#pragma once

#include <cmath>

#include <ddc/ddc.hpp>

#include <sll/mapping/curvilinear2d_to_cartesian.hpp>

template <class DimX, class DimY, class DimR, class DimP>
class CircularToCartesian : public Curvilinear2DToCartesian<DimX, DimY, DimR, DimP>
{
public:
    using cartesian_tag_x = DimX;
    using cartesian_tag_y = DimY;
    using circular_tag_r = DimR;
    using circular_tag_p = DimP;

public:
    CircularToCartesian() = default;

    CircularToCartesian(CircularToCartesian const& other) = default;

    CircularToCartesian(CircularToCartesian&& x) = default;

    ~CircularToCartesian() = default;

    CircularToCartesian& operator=(CircularToCartesian const& x) = default;

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
};
