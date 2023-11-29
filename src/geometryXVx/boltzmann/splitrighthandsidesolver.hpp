// SPDX-License-Identifier: MIT

#pragma once

#include <utility>
#include <vector>

#include <geometry.hpp>
#include <irighthandside.hpp>

#include "iboltzmannsolver.hpp"

/**
 * @brief A class that solves a Boltzmann equation using Strang's splitting.
 *
 * The solver splits the Boltzmann equation and separates the advective 
 * part from the source part. The sources refers to any operator that 
 * appears on the right-hand-side of Boltzmann's equation (typically, every 
 * operator except the advections). The splitting involves solving all the
 * source terms on a dt/2 timestep, then solving the advections on a dt
 * timestep using a Vlasov solver, then solving the sources again on dt/2
 * in reverse order. 
 */
class SplitRightHandSideSolver : public IBoltzmannSolver
{
    /** Member solver for the Vlasov equation. */
    IBoltzmannSolver const& m_boltzmann_solver;

    /** Member vector containing the source terms. */
    std::vector<std::reference_wrapper<IRightHandSide const>> m_rhs;

public:
    /**
     * @brief Creates an instance of the split boltzmann solver class.
     * @param[in] vlasov_solver A solver for the associated Vlasov equation 
     *                          (the boltzmann equation with no sources).
     * @param[in] rhs A vector containing all of the source terms of the 
     *                          considered Boltzmann equation.
     */
    SplitRightHandSideSolver(
            IBoltzmannSolver const& vlasov_solver,
            std::vector<std::reference_wrapper<IRightHandSide const>> rhs);

    ~SplitRightHandSideSolver() override = default;
    /**
     * @brief Solves a Boltzmann equation on a timestep dt.
     * @param[in, out] allfdistribu On input: the initial value of the distribution function.
     *                              On output: the value of the distribution function after solving 
     *                              the Boltzmann equation.
     * @param[in] electric_field The electric field computed at all spatial positions. 
     * @param[in] dt The timestep. 
     * @return The distribution function after solving the Boltzmann equation.
     */
    device_t<DSpanSpXVx> operator()(
            device_t<DSpanSpXVx> allfdistribu,
            DViewX electric_field,
            double dt) const override;
};
