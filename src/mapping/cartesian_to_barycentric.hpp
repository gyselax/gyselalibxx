// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "mapping_tools.hpp"

/**
 * @brief A class to convert cartesian coordinates to barycentric coordinates
 * on a triangle.
 *
 * Tags are used to identify the corners of the triangle. This ensures that
 * there are different types for coordinate systems related to different triangles.
 *
 * @tparam X The tag of the x dimension of the cartesian coordinates.
 * @tparam Y The tag of the y dimension of the cartesian coordinates.
 * @tparam Corner1Tag A tag identifying the first corner of the triangle.
 * @tparam Corner2Tag A tag identifying the second corner of the triangle.
 * @tparam Corner3Tag A tag identifying the third corner of the triangle.
 */
template <class X, class Y, class Corner1Tag, class Corner2Tag, class Corner3Tag>
class CartesianToBarycentric
{
public:
    /// The type of a coordinate in the barycentric coordinate system.
    using CoordResult = ddc::Coordinate<Corner1Tag, Corner2Tag, Corner3Tag>;

    /// The type of a coordinate in the cartesian coordinate system.
    using CoordArg = ddc::Coordinate<X, Y>;

private:
    using CartesianCoord = CoordArg;
    using BarycentricCoord = CoordResult;

private:
    CartesianCoord m_corner1;
    CartesianCoord m_corner2;
    CartesianCoord m_corner3;

public:
    /**
     * @brief Construct the operator which converts between the coordinate systems.
     *
     * @param[in] corner1 The coordinates of the first corner of the triangle.
     * @param[in] corner2 The coordinates of the second corner of the triangle.
     * @param[in] corner3 The coordinates of the third corner of the triangle.
     */
    CartesianToBarycentric(
            CartesianCoord const& corner1,
            CartesianCoord const& corner2,
            CartesianCoord const& corner3)
        : m_corner1(corner1)
        , m_corner2(corner2)
        , m_corner3(corner3)
    {
    }

    /**
     * @brief A copy operator for the mapping operator.
     * @param other The object to be copied.
     */
    CartesianToBarycentric(CartesianToBarycentric const& other) = default;

    /**
     * @brief A r-value copy operator for the mapping operator.
     * @param x The object to be consumed.
     */
    CartesianToBarycentric(CartesianToBarycentric&& x) = default;

    /// @brief The destructor of the mapping operator.
    ~CartesianToBarycentric() = default;

    /**
     * @brief A copy operator for the mapping operator.
     * @param x The object to be copied.
     * @return A reference to this class instance.
     */
    CartesianToBarycentric& operator=(CartesianToBarycentric const& x) = default;

    /**
     * @brief A r-value copy operator for the mapping operator.
     * @param x The object to be consumed.
     * @return A reference to this class instance.
     */
    CartesianToBarycentric& operator=(CartesianToBarycentric&& x) = default;

    /**
     * @brief The operator to get the equivalent barycentric coordinate of the cartesian coordinate.
     *
     * @param[in] pos The known cartesian coordinate.
     *
     * @return The equivalent barycentric coordinate.
     */
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
};
