// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

#include <iboltzmannsolver.hpp>
#include <ipoissonsolver.hpp>

#include "predcorr.hpp"

PredCorr::PredCorr(
        IBoltzmannSolver const& boltzmann_solver,
        IPoissonSolver const& poisson_solver,
        double const dt)
    : m_boltzmann_solver(boltzmann_solver)
    , m_poisson_solver(poisson_solver)
    , m_dt(dt)
{
}

DSpanSpXVx PredCorr::operator()(DSpanSpXVx const allfdistribu, int const steps) const
{
    // electrostatic potential and electric field (depending only on x)
    DFieldX electrostatic_potential(allfdistribu.domain<IDimX>());
    DFieldX electric_field(allfdistribu.domain<IDimX>());

    // a 2D chunck of the same size as fdistribu
    DFieldSpXVx allfdistribu_half_t(allfdistribu.domain());

    m_poisson_solver(electrostatic_potential, electric_field, allfdistribu);

    int iter = 0;
    for (; iter < steps; ++iter) {
        double const iter_time = iter * m_dt;

        // computation of the electrostatic potential at time tn and
        // the associated electric field
        m_poisson_solver(electrostatic_potential, electric_field, allfdistribu);

        PdiEvent("iteration")
                .with("iter", iter)
                .and_with("time_saved", iter_time)
                .and_with("fdistribu", allfdistribu)
                .and_with("electrostatic_potential", electrostatic_potential);

        // copy fdistribu
        deepcopy(allfdistribu_half_t, allfdistribu);

        // predictor
        m_boltzmann_solver(allfdistribu_half_t, electric_field, m_dt / 2);

        // computation of the electrostatic potential at time tn+1/2
        // and the associated electric field
        m_poisson_solver(electrostatic_potential, electric_field, allfdistribu_half_t);
        // correction on a dt
        m_boltzmann_solver(allfdistribu, electric_field, m_dt);
    }

    double const final_time = iter * m_dt;
    m_poisson_solver(electrostatic_potential, electric_field, allfdistribu);
    PdiEvent("last_iteration")
            .with("iter", iter)
            .and_with("time_saved", final_time)
            .and_with("fdistribu", allfdistribu)
            .and_with("electrostatic_potential", electrostatic_potential);

    return allfdistribu;
}
