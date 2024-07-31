// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "maxwellianequilibrium.hpp"

MaxwellianEquilibrium::MaxwellianEquilibrium(
        host_t<DFieldMemSp> density_eq,
        host_t<DFieldMemSp> temperature_eq,
        host_t<DFieldMemSp> mean_velocity_eq)
    : m_density_eq(std::move(density_eq))
    , m_temperature_eq(std::move(temperature_eq))
    , m_mean_velocity_eq(std::move(mean_velocity_eq))
{
}

DFieldSpVx MaxwellianEquilibrium::operator()(DFieldSpVx const allfequilibrium) const
{
    IdxRangeVx const gridvx = get_idx_range<GridVx>(allfequilibrium);
    IdxRangeSp const gridsp = get_idx_range<Species>(allfequilibrium);

    // Initialization of the maxwellian
    DFieldMemVx maxwellian_alloc(gridvx);
    ddc::ChunkSpan maxwellian = get_field(maxwellian_alloc);
    ddc::for_each(gridsp, [&](IdxSp const isp) {
        compute_maxwellian(
                maxwellian,
                m_density_eq(isp),
                m_temperature_eq(isp),
                m_mean_velocity_eq(isp));

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                gridvx,
                KOKKOS_LAMBDA(IdxVx const ivx) { allfequilibrium(isp, ivx) = maxwellian(ivx); });
    });
    return allfequilibrium;
}


MaxwellianEquilibrium MaxwellianEquilibrium::init_from_input(
        IdxRangeSp dom_kinsp,
        PC_tree_t const& yaml_input_file)
{
    host_t<DFieldMemSp> density_eq(dom_kinsp);
    host_t<DFieldMemSp> temperature_eq(dom_kinsp);
    host_t<DFieldMemSp> mean_velocity_eq(dom_kinsp);

    for (IdxSp const isp : dom_kinsp) {
        PC_tree_t const conf_isp = PCpp_get(yaml_input_file, ".SpeciesInfo[%d]", isp.uid());

        density_eq(isp) = PCpp_double(conf_isp, ".density_eq");
        temperature_eq(isp) = PCpp_double(conf_isp, ".temperature_eq");
        mean_velocity_eq(isp) = PCpp_double(conf_isp, ".mean_velocity_eq");
    }

    return MaxwellianEquilibrium(
            std::move(density_eq),
            std::move(temperature_eq),
            std::move(mean_velocity_eq));
}


void MaxwellianEquilibrium::compute_maxwellian(
        DFieldVx const fMaxwellian,
        double const density,
        double const temperature,
        double const mean_velocity)
{
    double const inv_sqrt_2piT = 1. / Kokkos::sqrt(2. * M_PI * temperature);
    IdxRangeVx const gridvx = get_idx_range(fMaxwellian);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridvx,
            KOKKOS_LAMBDA(IdxVx const ivx) {
                CoordVx const vx = ddc::coordinate(ivx);
                fMaxwellian(ivx) = density * inv_sqrt_2piT
                                   * Kokkos::exp(
                                           -(vx - mean_velocity) * (vx - mean_velocity)
                                           / (2. * temperature));
            });
}
