// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

/**
 * @brief An operator which solves the Poisson equation.
 *
 * An operator which solves the Poisson equation:
 * @f$ - \frac{d^2 \phi}{dx^2} = \rho @f$
 */
class IPoissonSolver
{
public:
    virtual ~IPoissonSolver() = default;

    /**
     * The operator which solves the equation using the method described by the class.
     *
     * @param[out] electrostatic_potential_device The electrostatic potential, the result of the poisson solver.
     * @param[out] electric_field_device The electric field, the derivative of the electrostatic potential.
     * @param[in] allfdistribu_device The distribution function.
     */
    virtual void operator()(
            device_t<DSpanX> electrostatic_potential_device,
            device_t<DSpanX> electric_field_device,
            device_t<DViewSpXVx> allfdistribu_device) const = 0;
};
