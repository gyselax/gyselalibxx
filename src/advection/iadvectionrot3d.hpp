// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

/**
 * @brief A class which provides an advection operator in a exact splitting for a shifted rotation.
 *
 * An abstract class which implements a function realizing the exact splitting for a shifted rotation.
 */
template <class Geometry, class GridVx, class GridVy, class GridVz>
class IAdvectionVelocityRot3D
{
public:
    virtual ~IAdvectionVelocityRot3D() = default;

    /**
     * @brief operates a transport of the distribution function.
     *
     * @param[in, out] allfdistribu Reference to an array containing the value of the distribution function.
     * @param[in] magnetic_field The magnetic field.
     * @param[in] mean_velocity_x The velocity frame shift in vx in the exact splitting.
     * @param[in] mean_velocity_y The velocity frame shift in vy in the exact splitting.
     * @param[in] dt Time step.
     *
     * @return A reference to an array containing the value of distribution the function at the updated time t+dt.
     */
    virtual DField<typename Geometry::IdxRangeFdistribu> operator()(
            DField<typename Geometry::IdxRangeFdistribu> allfdistribu,
            DConstField<typename Geometry::IdxRangeSpatial> magnetic_field_x,
            DConstField<typename Geometry::IdxRangeSpatial> magnetic_field_y,
            DConstField<typename Geometry::IdxRangeSpatial> magnetic_field_z,
            DConstField<typename Geometry::IdxRangeSpatial> mean_velocity_x,
            DConstField<typename Geometry::IdxRangeSpatial> mean_velocity_y,
            DConstField<typename Geometry::IdxRangeSpatial> mean_velocity_z,
            double dt) const = 0;
};
