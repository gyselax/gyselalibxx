// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

/**
 * @brief An operator which solves the Quasi-Neutrality equation.
 *
 * An operator which solves the Quasi-Neutrality equation:
 * @f$ - \frac{d^2 \phi}{dx^2} = \rho @f$
 */
class IQNSolver
{
public:
    virtual ~IQNSolver() = default;

    /**
     * The operator which solves the equation using the method described by the class.
     *
     * @param[out] electrostatic_potential The electrostatic potential, the result of the poisson solver.
     * @param[out] electric_field The electric field, the derivative of the electrostatic potential.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(
            DFieldX electrostatic_potential,
            DFieldX electric_field,
            DConstFieldSpXVx allfdistribu) const = 0;
};
