// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

#include <iboltzmannsolver.hpp>
#include <ipoissonsolver.hpp>

#include "predcorr.hpp"

PredCorr::PredCorr(IBoltzmannSolver const& boltzmann_solver, IPoissonSolver const& poisson_solver)
    : m_boltzmann_solver(boltzmann_solver)
    , m_poisson_solver(poisson_solver)
{
}

device_t<DSpanSpXVx> PredCorr::operator()(
        device_t<DSpanSpXVx> const allfdistribu_device,
        double const time_start,
        double const dt,
        int const steps) const
{
    auto allfdistribu_alloc = ddc::create_mirror_view_and_copy(allfdistribu_device);
    ddc::ChunkSpan allfdistribu = allfdistribu_alloc.span_view();

    // electrostatic potential and electric field (depending only on x)
    DFieldX electrostatic_potential(allfdistribu.domain<IDimX>());
    DFieldX electric_field(allfdistribu.domain<IDimX>());
    device_t<DFieldX> electric_field_device(allfdistribu.domain<IDimX>());

    // a 2D chunk of the same size as fdistribu
    DFieldSpXVx allfdistribu_half_t(allfdistribu.domain());
    device_t<DFieldSpXVx> allfdistribu_half_t_device(allfdistribu_device.domain());

    m_poisson_solver(electrostatic_potential, electric_field, allfdistribu);

    int iter = 0;
    for (; iter < steps; ++iter) {
        double const iter_time = time_start + iter * dt;

        // computation of the electrostatic potential at time tn and
        // the associated electric field
        ddc::deepcopy(allfdistribu, allfdistribu_device);
        m_poisson_solver(electrostatic_potential, electric_field, allfdistribu);

        ddc::PdiEvent("iteration")
                .with("iter", iter)
                .and_with("time_saved", iter_time)
                .and_with("fdistribu", allfdistribu)
                .and_with("electrostatic_potential", electrostatic_potential);

        // copy fdistribu
        ddc::deepcopy(allfdistribu_half_t_device, allfdistribu_device);
        ddc::deepcopy(electric_field_device, electric_field);

        // predictor
        m_boltzmann_solver(allfdistribu_half_t_device, electric_field_device, dt / 2);

        // computation of the electrostatic potential at time tn+1/2
        // and the associated electric field
        ddc::deepcopy(allfdistribu_half_t, allfdistribu_half_t_device);
        m_poisson_solver(electrostatic_potential, electric_field, allfdistribu_half_t);
        // correction on a dt
        ddc::deepcopy(electric_field_device, electric_field);
        m_boltzmann_solver(allfdistribu_device, electric_field_device, dt);
    }

    double const final_time = time_start + iter * dt;
    ddc::deepcopy(allfdistribu, allfdistribu_device);
    m_poisson_solver(electrostatic_potential, electric_field, allfdistribu);
    ddc::PdiEvent("last_iteration")
            .with("iter", iter)
            .and_with("time_saved", final_time)
            .and_with("fdistribu", allfdistribu)
            .and_with("electrostatic_potential", electrostatic_potential);

    return allfdistribu_device;
}
