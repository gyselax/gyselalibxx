// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "iqnsolver.hpp"
#include "ivlasovsolver.hpp"
#include "predcorr.hpp"
#include "transpose.hpp"

PredCorr::PredCorr(IVlasovSolver const& vlasov_solver, IQNSolver const& poisson_solver)
    : m_vlasov_solver(vlasov_solver)
    , m_poisson_solver(poisson_solver)
{
}

DFieldSpVxVyXY PredCorr::operator()(
        DFieldSpVxVyXY const allfdistribu_v2D_split,
        double const dt,
        int const steps) const
{
    IdxRangeSpXYVxVy idx_range_v2D_split_output_layout(get_idx_range(allfdistribu_v2D_split));
    DFieldMemSpXYVxVy allfdistribu_v2D_split_output_layout(idx_range_v2D_split_output_layout);
    auto allfdistribu_host_alloc
            = ddc::create_mirror_view(get_field(allfdistribu_v2D_split_output_layout));
    host_t<DFieldSpXYVxVy> allfdistribu_host = get_field(allfdistribu_host_alloc);

    // electrostatic potential and electric field (depending only on x)
    DFieldMemXY electrostatic_potential(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY electric_field_x(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY electric_field_y(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));

    host_t<DFieldMemXY> electrostatic_potential_host(
            get_idx_range<GridX, GridY>(allfdistribu_v2D_split));

    // a 2D memory block of the same size as fdistribu
    DFieldMemSpVxVyXY allfdistribu_half_t(get_idx_range(allfdistribu_v2D_split));

    m_poisson_solver(
            get_field(electrostatic_potential),
            get_field(electric_field_x),
            get_field(electric_field_y),
            get_const_field(allfdistribu_v2D_split));

    int iter = 0;
    for (; iter < steps; ++iter) {
        double const iter_time = iter * dt;

        // computation of the electrostatic potential at time tn and
        // the associated electric field
        m_poisson_solver(
                get_field(electrostatic_potential),
                get_field(electric_field_x),
                get_field(electric_field_y),
                get_const_field(allfdistribu_v2D_split));

        Kokkos::Profiling::pushRegion("PDIWrite");
        transpose_layout(
                Kokkos::DefaultExecutionSpace(),
                get_field(allfdistribu_v2D_split_output_layout),
                get_const_field(allfdistribu_v2D_split));
        // copies necessary to PDI
        ddc::parallel_deepcopy(
                allfdistribu_host,
                get_const_field(allfdistribu_v2D_split_output_layout));
        ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);
        ddc::PdiEvent("iteration")
                .with("iter", iter)
                .with("time_saved", iter_time)
                .with("fdistribu", allfdistribu_host)
                .with("electrostatic_potential", electrostatic_potential_host);
        // copy fdistribu
        ddc::parallel_deepcopy(allfdistribu_half_t, allfdistribu_v2D_split);
        Kokkos::Profiling::popRegion();

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
                get_field(allfdistribu_v2D_split),
                get_const_field(electric_field_x),
                get_const_field(electric_field_y),
                dt);
    }

    double const final_time = iter * dt;
    m_poisson_solver(
            get_field(electrostatic_potential),
            get_field(electric_field_x),
            get_field(electric_field_y),
            get_const_field(allfdistribu_v2D_split));

    Kokkos::Profiling::pushRegion("PDIWrite");
    transpose_layout(
            Kokkos::DefaultExecutionSpace(),
            get_field(allfdistribu_v2D_split_output_layout),
            get_const_field(allfdistribu_v2D_split));
    //copies necessary to PDI
    ddc::parallel_deepcopy(
            allfdistribu_host,
            get_const_field(allfdistribu_v2D_split_output_layout));
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);
    ddc::PdiEvent("last_iteration")
            .with("iter", iter)
            .with("time_saved", final_time)
            .with("fdistribu", allfdistribu_host)
            .with("electrostatic_potential", electrostatic_potential_host);
    Kokkos::Profiling::popRegion();

    return allfdistribu_v2D_split;
}
