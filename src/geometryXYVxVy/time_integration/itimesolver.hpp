// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for solving a Vlasov-Poisson system of equations.
 */
class ITimeSolver
{
public:
    virtual ~ITimeSolver() = default;

    /**
     * @brief Solves the Vlasov-Poisson system.
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving
     *                              the Vlasov-Poisson system a given number of iterations.
     * @param[in] dt The timestep.
     * @param[in] steps The number of iterations to be performed by the predictor-corrector.
     * @return The distribution function after solving the system.
     */
    virtual DFieldSpVxVyXY operator()(DFieldSpVxVyXY allfdistribu, double dt, int steps = 1)
            const = 0;
};
