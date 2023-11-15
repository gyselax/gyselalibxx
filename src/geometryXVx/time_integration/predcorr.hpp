// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "itimesolver.hpp"

class IPoissonSolver;
class IBoltzmannSolver;

/**
 * @brief A class that solves a Boltzmann-Poisson system of equations using a predictor-corrector scheme.
 *
 * A class that solves a Boltzmann-Poisson system with
 * a predictor corrector scheme. This scheme consists in 
 * estimating the electric potential after a time interval 
 * of a half-timestep. This potential is then used to compute
 * the value of the distribution function at time t+dt, where 
 * dt is the timestep.
 */
class PredCorr : public ITimeSolver
{
private:
    IBoltzmannSolver const& m_boltzmann_solver;

    IPoissonSolver const& m_poisson_solver;

public:
    /**
     * @brief Creates an instance of the predictor-corrector class.
     * @param[in] boltzmann_solver A solver for a Boltzmann equation.
     * @param[in] poisson_solver A solver for a Poisson equation.
     */
    PredCorr(IBoltzmannSolver const& boltzmann_solver, IPoissonSolver const& poisson_solver);

    ~PredCorr() override = default;

    /**
     * @brief Solves the Boltzmann-Poisson system.
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving 
     *                              the Boltzmann-Poisson system a given number of iterations.
     * @param[in] time_start The physical time at the start of the simulation.
     * @param[in] dt The timestep.
     * @param[in] steps The number of iterations to be performed by the predictor-corrector.
     * @return The distribution function after solving the system.
     */
    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double time_start, double dt, int steps = 1)
            const override;
};
