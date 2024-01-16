// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "ipoissonsolver.hpp"
/**
 * @brief Null operator.
 *
 * An operator which does not solve any equation.
 * It could be passed in operators which needs it, and used when solving poisson 
 * equation not relevant for the physics to study.
 */
class NullPoissonSolver : public IPoissonSolver
{
public:
    NullPoissonSolver() = default;

    ~NullPoissonSolver() override = default;
    /**
     * The operator which does not solves the equation.
     *
     * @param[out] electrostatic_potential_device The electrostatic potential, the result of the poisson solver.
     * @param[out] electric_field_device The electric field, the derivative of the electrostatic potential.
     * @param[in] allfdistribu_device The distribution function.
     */
    void operator()(
            device_t<DSpanX> const electrostatic_potential_device,
            device_t<DSpanX> const electric_field_device,
            device_t<DViewSpXVx> const allfdistribu_device) const override;
};
