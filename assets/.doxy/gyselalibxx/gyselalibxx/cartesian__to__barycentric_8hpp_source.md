

# File cartesian\_to\_barycentric.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**cartesian\_to\_barycentric.hpp**](cartesian__to__barycentric_8hpp.md)

[Go to the documentation of this file](cartesian__to__barycentric_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"

template <class Corner1Tag, class Corner2Tag, class Corner3Tag, class X, class Y>
class BarycentricToCartesian;

template <class X, class Y, class Corner1Tag, class Corner2Tag, class Corner3Tag>
class CartesianToBarycentric
{
public:
    using CoordResult = Coord<Corner1Tag, Corner2Tag, Corner3Tag>;

    using CoordArg = Coord<X, Y>;

private:
    using CartesianCoord = CoordArg;
    using BarycentricCoord = CoordResult;

private:
    CartesianCoord m_corner1;
    CartesianCoord m_corner2;
    CartesianCoord m_corner3;

public:
    CartesianToBarycentric(
            CartesianCoord const& corner1,
            CartesianCoord const& corner2,
            CartesianCoord const& corner3)
        : m_corner1(corner1)
        , m_corner2(corner2)
        , m_corner3(corner3)
    {
    }

    CartesianToBarycentric(CartesianToBarycentric const& other) = default;

    CartesianToBarycentric(CartesianToBarycentric&& x) = default;

    ~CartesianToBarycentric() = default;

    CartesianToBarycentric& operator=(CartesianToBarycentric const& x) = default;

    CartesianToBarycentric& operator=(CartesianToBarycentric&& x) = default;

    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& pos) const
    {
        const double x = ddc::get<X>(pos);
        const double y = ddc::get<Y>(pos);
        const double x1 = ddc::get<X>(m_corner1);
        const double x2 = ddc::get<X>(m_corner2);
        const double x3 = ddc::get<X>(m_corner3);
        const double y1 = ddc::get<Y>(m_corner1);
        const double y2 = ddc::get<Y>(m_corner2);
        const double y3 = ddc::get<Y>(m_corner3);

        const double div = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
        const double lam1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / div;
        const double lam2 = ((y - y3) * (x1 - x3) + (x3 - x) * (y1 - y3)) / div;
        const double lam3 = 1 - lam1 - lam2;
        return BarycentricCoord(lam1, lam2, lam3);
    }

    BarycentricToCartesian<Corner1Tag, Corner2Tag, Corner3Tag, X, Y> get_inverse_mapping() const
    {
        return BarycentricToCartesian<
                Corner1Tag,
                Corner2Tag,
                Corner3Tag,
                X,
                Y>(m_corner1, m_corner2, m_corner3);
    }
};
```


