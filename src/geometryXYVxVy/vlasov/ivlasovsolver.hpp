// SPDX-License-Identifier: MIT

#pragma once

#include "geometry_xyvxvy.hpp"

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
     * @param[in] efield The electric field computed at all spatial positions.
     * @param[in] dt The timestep. 
     *
     * @return The distribution function after solving the Vlasov equation.
     */
    virtual DFieldSpVxVyXY operator()(
            DFieldSpVxVyXY allfdistribu,
            DVectorConstFieldXY efield,
            double dt) const = 0;
};
