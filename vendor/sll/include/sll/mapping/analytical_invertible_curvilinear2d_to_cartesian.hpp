#pragma once

#include <array>

#include <ddc/ddc.hpp>

#include <sll/mapping/curvilinear2d_to_cartesian.hpp>

/**
 * @brief A class for describing analytical invertible curvilinear 2D mappings from the logical domain to the physical domain.
 *
 * We define a class of mapping not costly invertible (because analyticaly invertible).
 *
 * @see Curvilinear2DToCartesian
 */
template <class DimX, class DimY, class DimR, class DimP>
class AnalyticalInvertibleCurvilinear2DToCartesian
    : public Curvilinear2DToCartesian<DimX, DimY, DimR, DimP>
{
public:
    virtual ~AnalyticalInvertibleCurvilinear2DToCartesian() {};

    virtual ddc::Coordinate<DimX, DimY> operator()(
            ddc::Coordinate<DimR, DimP> const& coord) const = 0;

    /**
     * @brief Compute the logical coordinates from the physical coordinates.
     *
     * This class defined analytical invertible mappings which is not always the
     * case for a general mapping (see Curvilinear2DToCartesian::operator()).
     *
     * @param[in] coord
     * 			The coordinates in the physical domain.
     *
     * @return The coordinates in the logical domain.
     *
     * @see Curvilinear2DToCartesian::operator()
     */
    virtual ddc::Coordinate<DimR, DimP> operator()(
            ddc::Coordinate<DimX, DimY> const& coord) const = 0;
};
