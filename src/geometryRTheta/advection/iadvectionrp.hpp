// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

/**
 * @brief Define the base class of 2D advection operators in polar domain.
 */
class IAdvectionRP
{
public:
    virtual ~IAdvectionRP() = default;

    /**
     * @brief Advect a function along the advection field given on dt
     * with a given advection field along XY.
     *
     * @param[in, out] allfdistribu
     *      The function to be advected.
     * @param[in] advection_field
     *      The advection field along the physical domain axes, XY.
     * @param[in] dt
     *      The time step.
     *
     * @return A ChunkSpan to the advected function (allfdistribu).
     */
    virtual DSpanRP operator()(
            DSpanRP allfdistribu,
            VectorDViewRP<RDimX, RDimY> advection_field,
            double const dt) const = 0;

    /**
     * @brief Advect a function along the advection field given on dt
     * with a given advection field along RP.
     *
     * @param[in, out] allfdistribu
     *      The function to be advected.
     * @param[in] advection_field
     *      The advection field along the logical domain axes, RP.
     * @param[in] advection_field_xy_center
     *      The advection field along the physical domain axes, XY
     *      at the center point.
     * @param[in] dt
     *      The time step.
     *
     * @return A ChunkSpan to the advected function (allfdistribu).
     */
    virtual DSpanRP operator()(
            DSpanRP allfdistribu,
            VectorDViewRP<RDimR, RDimP> advection_field,
            CoordXY const& advection_field_xy_center,
            double const dt) const = 0;
};
