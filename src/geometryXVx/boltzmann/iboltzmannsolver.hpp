// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for solving a Boltzmann equation.
 */
class IBoltzmannSolver
{
public:
    virtual ~IBoltzmannSolver() = default;

    /**
     * @brief Operator for solving the Boltzmann equation on one timestep.
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving 
     *                              the Boltzmann equation on one timestep.
     * @param[in] efield The electric field computed at every spatial position.
     * @param[in] dt The timestep.
     * @return The distribution function after solving the Boltzmann equation.
     */
    virtual DFieldSpXVx operator()(DFieldSpXVx allfdistribu, DConstFieldX efield, double dt)
            const = 0;
};
