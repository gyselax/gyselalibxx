// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

#include <iboltzmannsolver.hpp>
#include <iqnsolver.hpp>

#include "predcorr.hpp"

PredCorr::PredCorr(IBoltzmannSolver const& boltzmann_solver, IQNSolver const& poisson_solver)
    : m_boltzmann_solver(boltzmann_solver)
    , m_poisson_solver(poisson_solver)
{
}

DSpanSpXVx PredCorr::operator()(
        DSpanSpXVx const allfdistribu,
        double const time_start,
        double const dt,
        int const steps) const
{
    auto allfdistribu_alloc = ddc::create_mirror_view_and_copy(allfdistribu);
    ddc::ChunkSpan allfdistribu_host = allfdistribu_alloc.span_view();

    // electrostatic potential and electric field (depending only on x)
    host_t<DFieldX> electrostatic_potential_host(allfdistribu.domain<IDimX>());
    DFieldX electrostatic_potential(allfdistribu.domain<IDimX>());

    DFieldX electric_field(allfdistribu.domain<IDimX>());

    // a 2D chunk of the same size as fdistribu
    host_t<DFieldSpXVx> allfdistribu_half_t_host(allfdistribu.domain());
    DFieldSpXVx allfdistribu_half_t(allfdistribu.domain());

    m_poisson_solver(electrostatic_potential, electric_field, allfdistribu);

    int iter = 0;
    for (; iter < steps; ++iter) {
        double const iter_time = time_start + iter * dt;

        // computation of the electrostatic potential at time tn and
        // the associated electric field
        m_poisson_solver(electrostatic_potential, electric_field, allfdistribu);
        // copies necessary to PDI
        ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
        ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);
        ddc::PdiEvent("iteration")
                .with("iter", iter)
                .and_with("time_saved", iter_time)
                .and_with("fdistribu", allfdistribu_host)
                .and_with("electrostatic_potential", electrostatic_potential_host);

        // copy fdistribu
        ddc::parallel_deepcopy(allfdistribu_half_t, allfdistribu);

        // predictor
        m_boltzmann_solver(allfdistribu_half_t, electric_field, dt / 2);

        // computation of the electrostatic potential at time tn+1/2
        // and the associated electric field
        m_poisson_solver(electrostatic_potential, electric_field, allfdistribu_half_t);
        // correction on a dt
        m_boltzmann_solver(allfdistribu, electric_field, dt);
    }

    double const final_time = time_start + iter * dt;
    m_poisson_solver(electrostatic_potential, electric_field, allfdistribu);
    //copies necessary to PDI
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);
    ddc::PdiEvent("last_iteration")
            .with("iter", iter)
            .and_with("time_saved", final_time)
            .and_with("fdistribu", allfdistribu_host)
            .and_with("electrostatic_potential", electrostatic_potential_host);

    return allfdistribu;
}
