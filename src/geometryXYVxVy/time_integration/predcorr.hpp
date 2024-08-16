// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "itimesolver.hpp"

class IQNSolver;
class IVlasovSolver;

/**
 * @brief A class that solves a Vlasov-Poisson system of equations using a predictor-corrector scheme.
 *
 * A class that solves a Vlasov-Poisson system with
 * a predictor corrector scheme. This scheme consists in
 * estimating the electric potential after a time interval
 * of a half-timestep. This potential is then used to compute
 * the value of the distribution function at time t+dt, where
 * dt is the timestep.
 */
class PredCorr : public ITimeSolver
{
private:
    IVlasovSolver const& m_vlasov_solver;

    IQNSolver const& m_poisson_solver;

public:
    /**
     * @brief Creates an instance of the predictor-corrector class.
     * @param[in] vlasov_solver A solver for a Boltzmann equation.
     * @param[in] poisson_solver A solver for a Poisson equation.
     */
    PredCorr(IVlasovSolver const& vlasov_solver, IQNSolver const& poisson_solver);

    ~PredCorr() override = default;

    /**
     * @brief Solves the Vlasov-Poisson system.
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving
     *                              the Vlasov-Poisson system a given number of iterations.
     * @param[in] dt The timestep.
     * @param[in] steps The number of iterations to be performed by the predictor-corrector.
     * @return The distribution function after solving the system.
     */
    DFieldSpXYVxVy operator()(DFieldSpXYVxVy allfdistribu, double dt, int steps = 1) const override;
};
