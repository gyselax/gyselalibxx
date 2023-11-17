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
     * @brief Advect a function along the advection field given on dt.
     *
     * @param[in, out] allfdistribu
     *      The function to be advected.
     * @param[in] advection_field
     *      The advection field.
     * @param[in] dt
     *      The time step.
     *
     * @return A ChunkSpan to the advected function (allfdistribu).
     */
    virtual DSpanRP operator()(
            DSpanRP allfdistribu,
            VectorDViewRP<RDimX, RDimY> advection_field,
            double const dt)
            = 0;
};
