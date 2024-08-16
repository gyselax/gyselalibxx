// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
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
     * @param[out] electric_field_x The x-component of the electric field, the gradient of the electrostatic potential.
     * @param[out] electric_field_y The y-component of the electric field, the gradient of the electrostatic potential.
     * @param[in] allfdistribu The distribution function.
     */
    void operator()(
            DFieldXY electrostatic_potential,
            DFieldXY electric_field_x,
            DFieldXY electric_field_y,
            DConstFieldSpXYVxVy allfdistribu) const override;
};
