// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

/**
 * @brief An abstract class for solving a Boltzmann-Poisson system of equations.
 */
class ITimeSolver
{
public:
    virtual ~ITimeSolver() = default;

    /**
     * @brief Operator for solving the Boltzmann-Poisson system.
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving 
     *                              the Boltzmann-Poisson system a given number of iterations.
     * @param[in] time_start The physical time at the start of the simulation.
     * @param[in] dt The timestep.
     * @param[in] steps The number of iterations to be performed by the solver.
     * @return The distribution function after solving the system.
     */
    virtual DSpanSpXVx operator()(
            DSpanSpXVx allfdistribu,
            double time_start,
            double dt,
            int steps = 1) const = 0;
};
