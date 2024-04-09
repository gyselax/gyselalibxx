// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

/**
 * @brief An abstract class for solving a fluid model.
 */
class IFluidSolver
{
public:
    virtual ~IFluidSolver() = default;

    /**
     * @brief Operator for solving the fluid model on one timestep.
     * @param[in, out] fluid_moments On input : a span referencing the fluid species.
     *                               On output : a span referencing the fluid species updated
     *                                after solving the fluid model on one timestep.
     * @param[in] allfdistribu A constant view referencing the distribution function.
     * @param[in] efield A constant view referencing the electric field.
     * @param[in] dt The timestep.
     * @return a span referencing the fluid species after solving the fluid model on one timestep.
     */
    virtual DSpanSpMX operator()(
            DSpanSpMX fluid_moments,
            DViewSpXVx allfdistribu,
            DViewX efield,
            double dt) const = 0;
};
