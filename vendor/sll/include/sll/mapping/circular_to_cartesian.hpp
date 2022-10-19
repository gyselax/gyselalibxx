#pragma once

#include <cmath>

#include <ddc/ddc.hpp>

template <class DimX, class DimY, class DimR, class DimP>
class CircularToCartesian
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

    Coordinate<DimX, DimY> operator()(Coordinate<DimR, DimP> const& coord) const
    {
        const double r = get<DimR>(coord);
        const double p = get<DimP>(coord);
        const double x = r * std::cos(p);
        const double y = r * std::sin(p);
        return Coordinate<DimX, DimY>(x, y);
    }

    Coordinate<DimR, DimP> operator()(Coordinate<DimX, DimY> const& coord) const
    {
        const double x = get<DimX>(coord);
        const double y = get<DimY>(coord);
        const double r = std::sqrt(x * x + y * y);
        const double p = std::atan2(y, x);
        return Coordinate<DimR, DimP>(r, p);
    }
};
