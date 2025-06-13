// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"

/**
 * @brief An operator which solves the Quasi-Neutrality equation using a fast
 * Fourier transform.
 *
 * An operator which solves the Quasi-Neutrality equation:
 * @f$ - \frac{d^2 \phi}{dx^2} = \rho @f$
 * using a fast Fourier transform on a periodic domain.
 * This operator only works for equidistant points.
 */
class IQNSolver
{
public:
    virtual ~IQNSolver() = default;

    /**
     * The operator which solves the equation using the method described by the class.
     *
     * @param[out] electrostatic_potential The electrostatic potential, the result of the poisson solver.
     * @param[out] electric_field The electric field, the gradient of the electrostatic potential.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(
            DFieldXY electrostatic_potential,
            DVectorFieldXY electric_field,
            DConstFieldSpVxVyXY allfdistribu) const = 0;
};
