// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

#include <ipoissonsolver.hpp>
#include <ivlasovsolver.hpp>

#include "predcorr.hpp"

PredCorr::PredCorr(IVlasovSolver const& vlasov_solver, IPoissonSolver const& poisson_solver)
    : m_vlasov_solver(vlasov_solver)
    , m_poisson_solver(poisson_solver)
{
}

DSpanSpXYVxVy PredCorr::operator()(
        DSpanSpXYVxVy const allfdistribu,
        double const dt,
        int const steps) const
{
    // electrostatic potential and electric field (depending only on x)
    DFieldXY electrostatic_potential(allfdistribu.domain<IDimX, IDimY>());
    DFieldXY electric_field_x(allfdistribu.domain<IDimX, IDimY>());
    DFieldXY electric_field_y(allfdistribu.domain<IDimX, IDimY>());

    // a 2D chunck of the same size as fdistribu
    DFieldSpXYVxVy allfdistribu_half_t(allfdistribu.domain());

    m_poisson_solver(electrostatic_potential, electric_field_x, electric_field_y, allfdistribu);

    int iter = 0;
    for (; iter < steps; ++iter) {
        double const iter_time = iter * dt;

        // computation of the electrostatic potential at time tn and
        // the associated electric field
        m_poisson_solver(electrostatic_potential, electric_field_x, electric_field_y, allfdistribu);

        ddc::PdiEvent("iteration")
                .with("iter", iter)
                .and_with("time_saved", iter_time)
                .and_with("fdistribu", allfdistribu)
                .and_with("electrostatic_potential", electrostatic_potential);

        // copy fdistribu
        ddc::parallel_deepcopy(allfdistribu_half_t, allfdistribu);

        // predictor
        m_vlasov_solver(allfdistribu_half_t, electric_field_x, electric_field_y, dt / 2);

        // computation of the electrostatic potential at time tn+1/2
        // and the associated electric field
        m_poisson_solver(
                electrostatic_potential,
                electric_field_x,
                electric_field_y,
                allfdistribu_half_t);

        // correction on a dt
        m_vlasov_solver(allfdistribu, electric_field_x, electric_field_y, dt);
    }

    double const final_time = iter * dt;
    m_poisson_solver(electrostatic_potential, electric_field_x, electric_field_y, allfdistribu);
    ddc::PdiEvent("last_iteration")
            .with("iter", iter)
            .and_with("time_saved", final_time)
            .and_with("fdistribu", allfdistribu)
            .and_with("electrostatic_potential", electrostatic_potential);

    return allfdistribu;
}
