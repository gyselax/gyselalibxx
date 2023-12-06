// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

/**
 * @brief A class which provides an advection operator.
 *
 * An abstract class which implements a function that 
 * applies the transport along a physical space direction of the phase space.
 */
template <class Geometry, class DDimX>
class IAdvectionSpatial
{
public:
    virtual ~IAdvectionSpatial() = default;
    /**
     * @brief operates a transport of the distribution function.
     *
     * @param[in, out] allfdistribu  reference to an array containing the value of distribution the function.
     * @param[in] dt time step.
     *
     * @return A reference to an array containing the value of distribution the function at the updated time T+dt.
     */
    virtual device_t<ddc::ChunkSpan<double, typename Geometry::FdistribuDDom>> operator()(
            device_t<ddc::ChunkSpan<double, typename Geometry::FdistribuDDom>> allfdistribu,
            double dt) const = 0;
};
