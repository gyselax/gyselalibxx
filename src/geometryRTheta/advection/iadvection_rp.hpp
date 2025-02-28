// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"

/**
 * @brief Define the base class of 2D advection operators in polar index range.
 */
class IAdvectionRTheta
{
public:
    virtual ~IAdvectionRTheta() = default;

    /**
     * @brief Advect a function along the advection field given on dt
     * with a given advection field along XY.
     *
     * @param[in, out] allfdistribu
     *      The function to be advected.
     * @param[in] advection_field
     *      The advection field along the physical index range axes, XY.
     * @param[in] dt
     *      The time step.
     *
     * @return A Field of the advected function (allfdistribu).
     */
    virtual host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> allfdistribu,
            host_t<DConstVectorFieldRTheta<X, Y>> advection_field,
            double const dt) const = 0;

    /**
     * @brief Advect a function along the advection field given on dt
     * with a given advection field along RTheta.
     *
     * @param[in, out] allfdistribu
     *      The function to be advected.
     * @param[in] advection_field
     *      The advection field along the logical index range axes, RTheta.
     * @param[in] advection_field_xy_centre
     *      The advection field along the physical index range axes, XY
     *      at the centre point.
     * @param[in] dt
     *      The time step.
     *
     * @return A Field of the advected function (allfdistribu).
     */
    virtual host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> allfdistribu,
            host_t<DConstVectorFieldRTheta<R, Theta>> advection_field,
            CoordXY const& advection_field_xy_centre,
            double const dt) const = 0;
};
