// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "maxwellianequilibrium.hpp"

MaxwellianEquilibrium::MaxwellianEquilibrium(
        DFieldSp density_eq,
        DFieldSp temperature_eq,
        DFieldSp mean_velocity_eq)
    : m_density_eq(std::move(density_eq))
    , m_temperature_eq(std::move(temperature_eq))
    , m_mean_velocity_eq(std::move(mean_velocity_eq))
{
}

DSpanSpVxVy MaxwellianEquilibrium::operator()(DSpanSpVxVy const allfequilibrium) const
{
    IDomainSp const gridsp = allfequilibrium.domain<IDimSp>();
    IDomainVx const gridvx = allfequilibrium.domain<IDimVx>();
    IDomainVy const gridvy = allfequilibrium.domain<IDimVy>();
    IDomainVxVy const gridvxvy = allfequilibrium.domain<IDimVx, IDimVy>();

    // Initialization of the maxwellian
    DFieldVxVy maxwellian(gridvxvy);
    ddc::for_each(gridsp, [&](IndexSp const isp) {
        compute_maxwellian(
                maxwellian,
                m_density_eq(isp),
                m_temperature_eq(isp),
                m_mean_velocity_eq(isp));

        ddc::for_each(gridvy, [&](IndexVy const ivy) {
            ddc::for_each(gridvx, [&](IndexVx const ivx) {
                allfequilibrium(isp, ivx, ivy) = maxwellian(ivx, ivy);
            });
        });
    });
    return allfequilibrium;
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
    IDomainVx const gridvx = fMaxwellian.domain<IDimVx>();
    IDomainVy const gridvy = fMaxwellian.domain<IDimVy>();

    for (IndexVx const ivx : gridvx) {
        double const vx = ddc::coordinate(ivx);
        for (IndexVy const ivy : gridvy) {
            double const vy = ddc::coordinate(ivy);
            fMaxwellian(ivx, ivy)
                    = density * inv_2pi
                      * std::exp(
                              -((vx - std::sqrt(mean_velocity)) * (vx - std::sqrt(mean_velocity))
                                + (vy - std::sqrt(mean_velocity)) * (vy - std::sqrt(mean_velocity)))
                              / (2. * temperature));
        }
    }
}
