// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

/**
 * @brief A class which provides an advection operator.
 *
 * An abstract class which implements a function that 
 * applies the transport along a velocity direction of the phase space.
 */
template <class Geometry, class GridV, class DataType>
class IAdvectionVelocity
{
public:
    virtual ~IAdvectionVelocity() = default;

    /**
     * @brief operates a transport of the distribution function.
     *
     * @param[in, out] allfdistribu Reference to an array containing the value of the distribution function.
     * @param[in] electric_field The electric field which derives from electrostatic potential and is the advection speed.
     * @param[in] dt Time step.
     *
     * @return A reference to an array containing the value of distribution the function at the updated time t+dt.
     */
    virtual Field<DataType, typename Geometry::IdxRangeFdistribu> operator()(
            Field<DataType, typename Geometry::IdxRangeFdistribu> allfdistribu,
            ConstField<DataType, typename Geometry::IdxRangeSpatial> electric_field,
            DataType dt) const = 0;
};
