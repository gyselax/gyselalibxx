

# File barycentric\_to\_cartesian.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**barycentric\_to\_cartesian.hpp**](barycentric__to__cartesian_8hpp.md)

[Go to the documentation of this file](barycentric__to__cartesian_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"

template <class Corner1Tag, class Corner2Tag, class Corner3Tag, class X, class Y>
class BarycentricToCartesian
{
public:
    using CoordArg = Coord<Corner1Tag, Corner2Tag, Corner3Tag>;

    using CoordResult = Coord<X, Y>;

private:
    using CartesianCoord = Coord<X, Y>;

private:
    CartesianCoord m_corner1;
    CartesianCoord m_corner2;
    CartesianCoord m_corner3;

public:
    BarycentricToCartesian(
            CartesianCoord const& corner1,
            CartesianCoord const& corner2,
            CartesianCoord const& corner3)
        : m_corner1(corner1)
        , m_corner2(corner2)
        , m_corner3(corner3)
    {
    }

    BarycentricToCartesian(BarycentricToCartesian const& other) = default;

    BarycentricToCartesian(BarycentricToCartesian&& x) = default;

    ~BarycentricToCartesian() = default;

    BarycentricToCartesian& operator=(BarycentricToCartesian const& x) = default;

    BarycentricToCartesian& operator=(BarycentricToCartesian&& x) = default;

    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& pos) const
    {
        const double l1 = ddc::get<Corner1Tag>(pos);
        const double l2 = ddc::get<Corner2Tag>(pos);
        const double l3 = ddc::get<Corner3Tag>(pos);

        const double x1 = ddc::get<X>(m_corner1);
        const double x2 = ddc::get<X>(m_corner2);
        const double x3 = ddc::get<X>(m_corner3);
        const double y1 = ddc::get<Y>(m_corner1);
        const double y2 = ddc::get<Y>(m_corner2);
        const double y3 = ddc::get<Y>(m_corner3);

        const double x = x1 * l1 + x2 * l2 + x3 * l3;
        const double y = y1 * l1 + y2 * l2 + y3 * l3;

        return CartesianCoord(x, y);
    }
};
```


