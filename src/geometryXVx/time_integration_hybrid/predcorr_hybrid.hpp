// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "ikineticfluidcoupling.hpp"
#include "itimesolver_hybrid.hpp"

class IQNSolver;
class IBoltzmannSolver;
class IFluidTransportSolver;
class IKineticFluidCoupling;

/**
 * @brief A class that solves a Boltzmann-Poisson system of equations coupled to a fluid particles model using a predictor-corrector scheme.
 *
 * A class that solves a Boltzmann-Poisson system with
 * a predictor corrector scheme. This scheme consists in 
 * estimating the electric potential after a time interval 
 * of a half-timestep. This potential is then used to compute
 * the value of the distribution function at time t+dt, and the
 * value of the fluid moments for the fluid species.
 * The fluid moments value is computed first without the source term S_n,N(x) via the fluid_solver.
 * Then, a coupling between the distribution function and the fluid moments is computed with the kinetic_fluid_coupling function.
 * dt is the timestep of the simulation.
 */
class PredCorrHybrid : public ITimeSolverHybrid
{
private:
    IBoltzmannSolver const& m_boltzmann_solver;

    IFluidTransportSolver const& m_fluid_solver;

    IQNSolver const& m_poisson_solver;

    IKineticFluidCoupling const& m_kinetic_fluid_coupling;

public:
    /**
     * @brief Creates an instance of the predictor-corrector class.
     * @param[in] boltzmann_solver A solver for a Boltzmann equation.
     * @param[in] fluid_solver A solver for a fluid model.
     * @param[in] poisson_solver A solver for a Quasi-Neutrality equation.
     * @param[in] kinetic_fluid_coupling A solver of the neutral source term in both the Boltzmann and fluid equations.
     */
    PredCorrHybrid(
            IBoltzmannSolver const& boltzmann_solver,
            IFluidTransportSolver const& fluid_solver,
            IQNSolver const& poisson_solver,
            IKineticFluidCoupling const& kinetic_fluid_coupling);

    ~PredCorrHybrid() override = default;

    /**
     * @brief Solves the Boltzmann-Poisson-fluid system.
     * @param[in, out] allfdistribu On input: the initial value of the distribution function.
     *                              On output: the value of the distribution function after solving 
     *                              the Boltzmann-Poisson-fluid system a given number of iterations.
     ** @param[in, out] fluid_moments On input: a span referencing the fluid species.
     *                                On output: the state of the fluid species after solving 
     *                                the Boltzmann-Poisson-fluid system a given number of iterations.
     * @param[in] time_start The physical time at the start of the simulation.
     * @param[in] dt The timestep.
     * @param[in] steps The number of iterations to be performed by the predictor-corrector.
     * @return The distribution function after solving the system.
     */
    DSpanSpXVx operator()(
            DSpanSpXVx allfdistribu,
            DSpanSpMX fluid_moments,
            double time_start,
            double dt,
            int steps = 1) const override;
};
