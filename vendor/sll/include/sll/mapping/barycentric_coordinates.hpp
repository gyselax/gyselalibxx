#pragma once

#include <ddc/ddc.hpp>

template <class DimX, class DimY, class Corner1Tag, class Corner2Tag, class Corner3Tag>
class CartesianToBarycentricCoordinates
{
public:
    using BarycentricCoord = ddc::Coordinate<Corner1Tag, Corner2Tag, Corner3Tag>;

    using CartesianCoord = ddc::Coordinate<DimX, DimY>;

private:
    CartesianCoord m_corner1;
    CartesianCoord m_corner2;
    CartesianCoord m_corner3;

public:
    CartesianToBarycentricCoordinates(
            CartesianCoord const& corner1,
            CartesianCoord const& corner2,
            CartesianCoord const& corner3)
        : m_corner1(corner1)
        , m_corner2(corner2)
        , m_corner3(corner3)
    {
    }

    CartesianToBarycentricCoordinates(CartesianToBarycentricCoordinates const& other) = default;

    CartesianToBarycentricCoordinates(CartesianToBarycentricCoordinates&& x) = default;

    ~CartesianToBarycentricCoordinates() = default;

    CartesianToBarycentricCoordinates& operator=(CartesianToBarycentricCoordinates const& x)
            = default;

    CartesianToBarycentricCoordinates& operator=(CartesianToBarycentricCoordinates&& x) = default;

    BarycentricCoord operator()(CartesianCoord const& pos) const
    {
        const double x = ddc::get<DimX>(pos);
        const double y = ddc::get<DimY>(pos);
        const double x1 = ddc::get<DimX>(m_corner1);
        const double x2 = ddc::get<DimX>(m_corner2);
        const double x3 = ddc::get<DimX>(m_corner3);
        const double y1 = ddc::get<DimY>(m_corner1);
        const double y2 = ddc::get<DimY>(m_corner2);
        const double y3 = ddc::get<DimY>(m_corner3);

        const double div = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
        const double lam1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / div;
        const double lam2 = ((y - y3) * (x1 - x3) + (x3 - x) * (y1 - y3)) / div;
        const double lam3 = 1 - lam1 - lam2;
        return BarycentricCoord(lam1, lam2, lam3);
    }

    CartesianCoord operator()(BarycentricCoord const& pos) const
    {
        const double l1 = ddc::get<Corner1Tag>(pos);
        const double l2 = ddc::get<Corner2Tag>(pos);
        const double l3 = ddc::get<Corner3Tag>(pos);

        const double x1 = ddc::get<DimX>(m_corner1);
        const double x2 = ddc::get<DimX>(m_corner2);
        const double x3 = ddc::get<DimX>(m_corner3);
        const double y1 = ddc::get<DimY>(m_corner1);
        const double y2 = ddc::get<DimY>(m_corner2);
        const double y3 = ddc::get<DimY>(m_corner3);

        const double x = x1 * l1 + x2 * l2 + x3 * l3;
        const double y = y1 * l1 + y2 * l2 + y3 * l3;

        return CartesianCoord(x, y);
    }
};
