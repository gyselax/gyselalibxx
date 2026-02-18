// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "geometry_r_theta.hpp"

/**
 * @brief Base class for a Quasi-Neutrality solver.
 */
class IQNSolver
{
public:
    virtual ~IQNSolver() = default;

    /**
     * @brief Compute the electrical potential and
     * the electric field from the Quasi-Neutrality equation.
     *
     * @param[out] electrostatic_potential
     *      The solution of the Quasi-Neutrality equation.
     * @param[out] electric_field
     *      The electric field @f$E = -\nabla \phi@f$.
     * @param[in] density
     *      The rhs of the Quasi-Neutrality equation.
     */
    virtual void operator()(
            host_t<DFieldRTheta> electrostatic_potential,
            host_t<DVectorFieldRTheta<X, Y>> electric_field,
            host_t<DConstFieldRTheta> density) const = 0;
};
