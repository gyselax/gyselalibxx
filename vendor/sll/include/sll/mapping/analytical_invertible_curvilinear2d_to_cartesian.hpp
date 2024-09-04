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
template <class X, class Y, class R, class Theta>
class AnalyticalInvertibleCurvilinear2DToCartesian : public Curvilinear2DToCartesian<X, Y, R, Theta>
{
public:
    KOKKOS_FUNCTION virtual ~AnalyticalInvertibleCurvilinear2DToCartesian() {}

    KOKKOS_FUNCTION virtual ddc::Coordinate<X, Y> operator()(
            ddc::Coordinate<R, Theta> const& coord) const = 0;

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
    KOKKOS_FUNCTION virtual ddc::Coordinate<R, Theta> operator()(
            ddc::Coordinate<X, Y> const& coord) const = 0;
};
