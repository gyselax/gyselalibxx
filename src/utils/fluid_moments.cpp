// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <fluid_moments.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

FluidMoments::FluidMoments(Quadrature<IDimVx> integrate_v) : m_integrate_v(std::move(integrate_v))
{
}

/*
 * Computes the density of allfdistribu
*/
void FluidMoments::operator()(
        DSpanSpX const density,
        DViewSpXVx const allfdistribu,
        FluidMoments::MomentDensity)
{
    for_each(
            policies::parallel_host,
            get_domain<IDimSp, IDimX>(allfdistribu),
            [&](IndexSpX const ispx) { density(ispx) = m_integrate_v(allfdistribu[ispx]); });
}

/*
 * Computes the mean_velocity of allfdistribu, using its density
*/
void FluidMoments::operator()(
        DSpanSpX const mean_velocity,
        DViewSpXVx const allfdistribu,
        DViewSpX const density,
        FluidMoments::MomentVelocity)
{
    DFieldSpXVx integrand(allfdistribu.domain());
    for_each(policies::parallel_host, allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        CoordVx const coordv = coordinate(select<IDimVx>(ispxvx));
        integrand(ispxvx) = coordv * allfdistribu(ispxvx);
    });

    for_each(
            policies::parallel_host,
            get_domain<IDimSp, IDimX>(allfdistribu),
            [&](IndexSpX const ispx) {
                mean_velocity(ispx) = m_integrate_v(integrand[ispx]) / density(ispx);
            });
}

/*
 * Computes the temperature of allfdistribu, using its density and mean velocity
*/
void FluidMoments::operator()(
        DSpanSpX const temperature,
        DViewSpXVx const allfdistribu,
        DViewSpX const density,
        DViewSpX const mean_velocity,
        FluidMoments::MomentTemperature)
{
    DFieldSpXVx integrand(allfdistribu.domain());
    for_each(policies::parallel_host, allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        double const coeff
                = coordinate(select<IDimVx>(ispxvx)) - mean_velocity(select<IDimSp, IDimX>(ispxvx));
        integrand(ispxvx) = coeff * coeff * allfdistribu(ispxvx);
    });

    for_each(
            policies::parallel_host,
            get_domain<IDimSp, IDimX>(allfdistribu),
            [&](IndexSpX const ispx) {
                temperature(ispx) = m_integrate_v(integrand[ispx]) / density(ispx);
            });
}
