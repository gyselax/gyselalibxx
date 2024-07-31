// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

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
     * @param[in] allfdistribu
     *      The rhs of the Quasi-Neutrality equation.
     */
    virtual void operator()(
            DFieldRTheta electrostatic_potential,
            DVectorFieldRTheta<X, Y> electric_field,
            DConstFieldRTheta allfdistribu) const = 0;
};
