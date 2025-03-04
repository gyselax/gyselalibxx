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

DFieldSpVxVy MaxwellianEquilibrium::operator()(DFieldSpVxVy const allfequilibrium) const
{
    IdxRangeSp const gridsp = get_idx_range<Species>(allfequilibrium);
    IdxRangeVxVy const gridvxvy = get_idx_range<GridVx, GridVy>(allfequilibrium);

    // Initialisation of the maxwellian
    DFieldMemVxVy maxwellian_alloc(gridvxvy);
    DFieldVxVy maxwellian = get_field(maxwellian_alloc);
    ddc::for_each(gridsp, [&](IdxSp const isp) {
        compute_maxwellian(
                maxwellian,
                m_density_eq(isp),
                m_temperature_eq(isp),
                m_mean_velocity_eq(isp));

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                gridvxvy,
                KOKKOS_LAMBDA(IdxVxVy const ivxvy) {
                    allfequilibrium(isp, ivxvy) = maxwellian(ivxvy);
                });
    });
    return allfequilibrium;
}


MaxwellianEquilibrium MaxwellianEquilibrium::init_from_input(
        IdxRangeSp idx_range_kinsp,
        PC_tree_t const& yaml_input_file)
{
    host_t<DFieldMemSp> density_eq(idx_range_kinsp);
    host_t<DFieldMemSp> temperature_eq(idx_range_kinsp);
    host_t<DFieldMemSp> mean_velocity_eq(idx_range_kinsp);

    for (IdxSp const isp : idx_range_kinsp) {
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


/*
 Computing the Maxwellian function as
  fM(vx,vy) = n/(2*PI*T)*exp(-(vx**2+vy**2)/(2*T))
 with n the density and T the temperature and
*/
void MaxwellianEquilibrium::compute_maxwellian(
        DFieldVxVy const fMaxwellian,
        double const density,
        double const temperature,
        double const mean_velocity)
{
    double const inv_2pi = 1. / (2. * M_PI * temperature);
    IdxRangeVxVy const gridvxvy = get_idx_range<GridVx, GridVy>(fMaxwellian);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridvxvy,
            KOKKOS_LAMBDA(IdxVxVy const ivxvy) {
                double const vx = ddc::coordinate(ddc::select<GridVx>(ivxvy));
                double const vy = ddc::coordinate(ddc::select<GridVy>(ivxvy));
                fMaxwellian(ivxvy) = density * inv_2pi
                                     * Kokkos::exp(
                                             -((vx - mean_velocity) * (vx - mean_velocity)
                                               + (vy - mean_velocity) * (vy - mean_velocity))
                                             / (2. * temperature));
            });
}
