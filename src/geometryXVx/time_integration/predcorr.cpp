// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "iboltzmannsolver.hpp"
#include "iqnsolver.hpp"
#include "predcorr.hpp"

PredCorr::PredCorr(IBoltzmannSolver const& boltzmann_solver, IQNSolver const& poisson_solver)
    : m_boltzmann_solver(boltzmann_solver)
    , m_poisson_solver(poisson_solver)
{
}

DFieldSpXVx PredCorr::operator()(
        DFieldSpXVx const allfdistribu,
        double const time_start,
        double const dt,
        int const steps) const
{
    auto allfdistribu_alloc = ddc::create_mirror_view(allfdistribu);
    host_t<DFieldSpXVx> allfdistribu_host = get_field(allfdistribu_alloc);

    // electrostatic potential and electric field (depending only on x)
    host_t<DFieldMemX> electrostatic_potential_host(get_idx_range<GridX>(allfdistribu));
    DFieldMemX electrostatic_potential(get_idx_range<GridX>(allfdistribu));

    DFieldMemX electric_field(get_idx_range<GridX>(allfdistribu));

    // a 2D chunk of the same size as fdistribu
    host_t<DFieldMemSpXVx> allfdistribu_half_t_host(get_idx_range(allfdistribu));
    DFieldMemSpXVx allfdistribu_half_t(get_idx_range(allfdistribu));

    m_poisson_solver(
            get_field(electrostatic_potential),
            get_field(electric_field),
            get_const_field(allfdistribu));

    int iter = 0;
    for (; iter < steps; ++iter) {
        Kokkos::Profiling::pushRegion("Time step");
        double const iter_time = time_start + iter * dt;

        // computation of the electrostatic potential at time tn and
        // the associated electric field
        m_poisson_solver(
                get_field(electrostatic_potential),
                get_field(electric_field),
                get_const_field(allfdistribu));
        // copies necessary to PDI
        ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
        ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);
        Kokkos::Profiling::pushRegion("HDF5_Output");
        ddc::PdiEvent("iteration")
                .with("iter", iter)
                .with("time_saved", iter_time)
                .with("fdistribu", allfdistribu_host)
                .with("electrostatic_potential", electrostatic_potential_host);
        Kokkos::Profiling::popRegion();

        // copy fdistribu
        ddc::parallel_deepcopy(allfdistribu_half_t, allfdistribu);

        // predictor
        m_boltzmann_solver(get_field(allfdistribu_half_t), get_const_field(electric_field), dt / 2);

        // computation of the electrostatic potential at time tn+1/2
        // and the associated electric field
        m_poisson_solver(
                get_field(electrostatic_potential),
                get_field(electric_field),
                get_const_field(allfdistribu_half_t));
        // correction on a dt
        m_boltzmann_solver(allfdistribu, get_const_field(electric_field), dt);

        Kokkos::Profiling::popRegion();
    }

    double const final_time = time_start + iter * dt;
    m_poisson_solver(
            get_field(electrostatic_potential),
            get_field(electric_field),
            get_const_field(allfdistribu));
    //copies necessary to PDI
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);
    ddc::PdiEvent("last_iteration")
            .with("iter", iter)
            .with("time_saved", final_time)
            .with("fdistribu", allfdistribu_host)
            .with("electrostatic_potential", electrostatic_potential_host);

    return allfdistribu;
}
