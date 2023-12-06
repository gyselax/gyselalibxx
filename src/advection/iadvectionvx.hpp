// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>

/**
 * @brief A class which provides an advection operator.
 *
 * An abstract class which implements a function that 
 * applies the transport along a velocity direction of the phase space.
 */
template <class Geometry, class DDimV>
class IAdvectionVelocity
{
public:
    virtual ~IAdvectionVelocity() = default;

    /**
     * @brief operates a transport of the distribution function.
     *
     * @param[in, out] allfdistribu Reference to an array containing the value of the distribution function.
     * @param[in] electrostatic_potential The electrostatic potential which is the advection speed.
     * @param[in] dt Time step.
     *
     * @return A reference to an array containing the value of distribution the function at the updated time t+dt.
     */
    virtual device_t<ddc::ChunkSpan<double, typename Geometry::FdistribuDDom>> operator()(
            device_t<ddc::ChunkSpan<double, typename Geometry::FdistribuDDom>> allfdistribu,
            device_t<ddc::ChunkSpan<const double, typename Geometry::SpatialDDom>>
                    electrostatic_potential,
            double dt) const = 0;
};
