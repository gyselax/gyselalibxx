// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

#include <iboltzmannsolver.hpp>
#include <ifluidtransportsolver.hpp>
#include <iqnsolver.hpp>

#include "predcorr_hybrid.hpp"

PredCorrHybrid::PredCorrHybrid(
        IBoltzmannSolver const& boltzmann_solver,
        IFluidTransportSolver const& fluid_solver,
        IQNSolver const& poisson_solver,
        IKineticFluidCoupling const& kinetic_fluid_coupling)
    : m_boltzmann_solver(boltzmann_solver)
    , m_fluid_solver(fluid_solver)
    , m_poisson_solver(poisson_solver)
    , m_kinetic_fluid_coupling(kinetic_fluid_coupling)
{
}

DFieldSpXVx PredCorrHybrid::operator()(
        DFieldSpXVx const allfdistribu,
        DFieldSpMomX fluid_moments,
        double const time_start,
        double const dt,
        int const steps) const
{
    auto allfdistribu_alloc = ddc::create_mirror_view_and_copy(allfdistribu);
    ddc::ChunkSpan allfdistribu_host = get_field(allfdistribu_alloc);

    IdxRangeX const dom_x = get_idx_range<GridX>(allfdistribu);

    // electrostatic potential and electric field (depending only on x)
    host_t<DFieldMemX> electrostatic_potential_host(dom_x);
    DFieldMemX electrostatic_potential(dom_x);

    DFieldMemX electric_field(dom_x);

    host_t<DFieldMemSpMomX> fluid_moments_host(get_idx_range(fluid_moments));

    // a 2D chunk of the same size as fdistribu
    host_t<DFieldMemSpXVx> allfdistribu_half_t_host(get_idx_range(allfdistribu));
    DFieldMemSpXVx allfdistribu_half_t(get_idx_range(allfdistribu));

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
        ddc::parallel_deepcopy(fluid_moments_host, fluid_moments);
        ddc::PdiEvent("iteration")
                .with("iter", iter)
                .and_with("time_saved", iter_time)
                .and_with("fdistribu", allfdistribu_host)
                .and_with("fluid_moments", fluid_moments_host)
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
        m_fluid_solver(fluid_moments, allfdistribu, electric_field, dt);

        m_kinetic_fluid_coupling(allfdistribu, fluid_moments, dt);
    }

    double const final_time = time_start + iter * dt;
    m_poisson_solver(electrostatic_potential, electric_field, allfdistribu);

    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);
    ddc::parallel_deepcopy(fluid_moments_host, fluid_moments);
    ddc::PdiEvent("last_iteration")
            .with("iter", iter)
            .and_with("time_saved", final_time)
            .and_with("fdistribu", allfdistribu_host)
            .and_with("fluid_moments", fluid_moments_host)
            .and_with("electrostatic_potential", electrostatic_potential_host);

    return allfdistribu;
}
