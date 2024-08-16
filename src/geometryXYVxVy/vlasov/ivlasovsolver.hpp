// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for solving a Vlasov equation.
 */
class IVlasovSolver
{
public:
    virtual ~IVlasovSolver() = default;

    /**
     * @brief Solves a Vlasov equation on a timestep dt.
     *
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving 
     *                              the Vlasov equation.
     * @param[in] efield_x The electric field in the x direction computed at all spatial positions. 
     * @param[in] efield_y The electric field in the y direction computed at all spatial positions. 
     * @param[in] dt The timestep. 
     *
     * @return The distribution function after solving the Vlasov equation.
     */
    virtual DFieldSpXYVxVy operator()(
            DFieldSpXYVxVy allfdistribu,
            DConstFieldXY efield_x,
            DConstFieldXY efield_y,
            double dt) const = 0;
};
