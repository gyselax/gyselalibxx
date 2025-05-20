// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "ihybridsplitting.hpp"

class IHybridFieldSolver;
class IHybridVlasovSolver;

/**
 * @brief A class that solves a hybrid plasma model with kinetic ions and massless electrons
 *        using a Strang splitting.
 *
 * A class that solves a hybrid plasma model with kinetic ions and massless electrons using
 * a Strang splitting. This scheme consists in a modified mid-point rule with a specially
 * chosen mean-velocity. The nonlinear iteration is only for the magnetic field, pressure, 
 * and the mean-velocity, after which the distribution functions are updated to t+dt by the 
 * updated magnetic field, pressure, and the mean-velocity.
 */
class Hybridsplitting : public IHybridSplitting
{
private:
    IHybridVlasovSolver const& m_vlasov_solver;

    IHybridFieldSolver const& m_hybrid_solver;

    IMomentsCalculator const& m_moments_calculator;

public:
    /**
     * @brief Creates an instance of the Strang splitting of the hybrid model.
     * @param[in] vlasov_solver A solver for a Boltzmann equation.
     * @param[in] hybrid_solver A solver for a Poisson equation.
     * @param[in] moments_calculator A solver for calculating the velocity moments.
     */
    Hybridsplitting(IHybridVlasovSolver const& vlasov_solver, IHybridFieldSolver const& hybrid_solver,
                    IMomentsCalculator const& moments_calculator);

    ~Hybridsplitting() override = default;

    /**
     * @brief Solves the hybrid model with kinetic ions and massless electrons.
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving
     *                              the hybrid model a given number of iterations.
     * @param[in] mean_velocity_x_each The mean velocity in vx for each species.
     * @param[in] mean_velocity_y_each The mean velocity in vy for each species.
     * @param[in] mean_velocity_x The mean velocity in vx for all species.
     * @param[in] momentum_x The momentum in vx for all species.
     * @param[in] momentum_y The momentum in vy for all species.
     * @param[in, out] magnetic_field_z The magnetic field.
     * @param[in, out] pressure The pressure field.
     * @param[in] rho_each The charge density for each species.
     * @param[in] rho The charge density for all species.
     * @param[out] kinetic The kinetic density for all species.
     * @param[in] dt The timestep.
     * @param[in] steps The number of iterations to be performed by the predictor-corrector.
     * @return The distribution function after solving the system.
     */
    DFieldSpVxVyXY operator()(DFieldSpVxVyXY allfdistribu, 
                              DFieldSpXY mean_velocity_x_each, 
                              DFieldSpXY mean_velocity_y_each,
                              DFieldXY mean_velocity_x,
                              DFieldXY mean_velocity_y,
                              DFieldXY momentum_x,
                              DFieldXY momentum_y,
                              DFieldXY magnetic_field_z,
                              DFieldXY pressure,
                              DFieldSpXY rho_each, 
                              DFieldXY rho,
                              DFieldXY kinetic, 
                              double dt, int steps = 1) const override;

};
