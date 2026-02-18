// SPDX-License-Identifier: MIT

#pragma once

#include "geometry_xyvxvy.hpp"
#include "iqnsolver.hpp"

/**
 * @brief A Quasi-Neutrality solver which does nothing.
 * This can potentially be useful when debugging.
 */
class NullQNSolver : public IQNSolver
{
public:
    NullQNSolver() = default;

    ~NullQNSolver() override = default;

    /** @brief A QN Solver which does nothing
     *
     * @param[out] electrostatic_potential The electrostatic potential, the result of the poisson solver.
     * @param[out] electric_field The electric field, the gradient of the electrostatic potential.
     * @param[in] allfdistribu The distribution function.
     */
    void operator()(
            DFieldXY electrostatic_potential,
            DVectorFieldXY electric_field,
            DConstFieldSpVxVyXY allfdistribu) const override;
};
