// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "iqnsolver.hpp"
/**
 * @brief Null operator.
 *
 * An operator which does not solve any equation.
 * It could be passed in operators which needs it, and used when solving poisson 
 * equation not relevant for the physics to study.
 */
class NullQNSolver : public IQNSolver
{
public:
    NullQNSolver() = default;

    ~NullQNSolver() override = default;
    /**
     * The operator which does not solves the equation.
     *
     * @param[out] electrostatic_potential The electrostatic potential, the result of the poisson solver.
     * @param[out] electric_field The electric field, the derivative of the electrostatic potential.
     * @param[in] allfdistribu The distribution function.
     */
    void operator()(
            DFieldX const electrostatic_potential,
            DFieldX const electric_field,
            DConstFieldSpXVx const allfdistribu) const override;
};
