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

DSpanSpVxVy MaxwellianEquilibrium::operator()(DSpanSpVxVy const allfequilibrium) const
{
    IdxRangeSp const gridsp = allfequilibrium.domain<Species>();
    IDomainVxVy const gridvxvy = allfequilibrium.domain<IDimVx, IDimVy>();

    // Initialization of the maxwellian
    DFieldVxVy maxwellian_alloc(gridvxvy);
    DSpanVxVy maxwellian = maxwellian_alloc.span_view();
    ddc::for_each(gridsp, [&](IdxSp const isp) {
        compute_maxwellian(
                maxwellian,
                m_density_eq(isp),
                m_temperature_eq(isp),
                m_mean_velocity_eq(isp));

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                gridvxvy,
                KOKKOS_LAMBDA(IndexVxVy const ivxvy) {
                    allfequilibrium(isp, ivxvy) = maxwellian(ivxvy);
                });
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


/*
 Computing the Maxwellian function as
  fM(vx,vy) = n/(2*PI*T)*exp(-(vx**2+vy**2)/(2*T))
 with n the density and T the temperature and
*/
void MaxwellianEquilibrium::compute_maxwellian(
        DSpanVxVy const fMaxwellian,
        double const density,
        double const temperature,
        double const mean_velocity)
{
    double const inv_2pi = 1. / (2. * M_PI * temperature);
    IDomainVxVy const gridvxvy = fMaxwellian.domain<IDimVx, IDimVy>();

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridvxvy,
            KOKKOS_LAMBDA(IndexVxVy const ivxvy) {
                double const vx = ddc::coordinate(ddc::select<IDimVx>(ivxvy));
                double const vy = ddc::coordinate(ddc::select<IDimVy>(ivxvy));
                fMaxwellian(ivxvy) = density * inv_2pi
                                     * Kokkos::exp(
                                             -((vx - mean_velocity) * (vx - mean_velocity)
                                               + (vy - mean_velocity) * (vy - mean_velocity))
                                             / (2. * temperature));
            });
}
