// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "iqnsolver.hpp"
#include "ivlasovsolver.hpp"
#include "predcorr.hpp"

PredCorr::PredCorr(IVlasovSolver const& vlasov_solver, IQNSolver const& poisson_solver)
    : m_vlasov_solver(vlasov_solver)
    , m_poisson_solver(poisson_solver)
{
}

DFieldSpXYVxVy PredCorr::operator()(
        DFieldSpXYVxVy const allfdistribu,
        double const dt,
        int const steps) const
{
    auto allfdistribu_host_alloc = ddc::create_mirror_view_and_copy(allfdistribu);
    host_t<DFieldSpXYVxVy> allfdistribu_host = get_field(allfdistribu_host_alloc);

    // electrostatic potential and electric field (depending only on x)
    DFieldMemXY electrostatic_potential(get_idx_range<GridX, GridY>(allfdistribu));
    DFieldMemXY electric_field_x(get_idx_range<GridX, GridY>(allfdistribu));
    DFieldMemXY electric_field_y(get_idx_range<GridX, GridY>(allfdistribu));

    host_t<DFieldMemXY> electrostatic_potential_host(get_idx_range<GridX, GridY>(allfdistribu));

    // a 2D chunck of the same size as fdistribu
    DFieldMemSpXYVxVy allfdistribu_half_t(get_idx_range(allfdistribu));

    m_poisson_solver(
            get_field(electrostatic_potential),
            get_field(electric_field_x),
            get_field(electric_field_y),
            get_const_field(allfdistribu));

    int iter = 0;
    for (; iter < steps; ++iter) {
        double const iter_time = iter * dt;

        // computation of the electrostatic potential at time tn and
        // the associated electric field
        m_poisson_solver(
                get_field(electrostatic_potential),
                get_field(electric_field_x),
                get_field(electric_field_y),
                get_const_field(allfdistribu));
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
        m_vlasov_solver(
                get_field(allfdistribu_half_t),
                get_const_field(electric_field_x),
                get_const_field(electric_field_y),
                dt / 2);

        // computation of the electrostatic potential at time tn+1/2
        // and the associated electric field
        m_poisson_solver(
                get_field(electrostatic_potential),
                get_field(electric_field_x),
                get_field(electric_field_y),
                get_const_field(allfdistribu_half_t));

        // correction on a dt
        m_vlasov_solver(
                get_field(allfdistribu),
                get_const_field(electric_field_x),
                get_const_field(electric_field_y),
                dt);
    }

    double const final_time = iter * dt;
    m_poisson_solver(
            get_field(electrostatic_potential),
            get_field(electric_field_x),
            get_field(electric_field_y),
            get_const_field(allfdistribu));

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
