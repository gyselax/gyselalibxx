// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <fluid_moments.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

FluidMoments::FluidMoments(Quadrature<IDimVx> integrate_v) : m_integrate_v(std::move(integrate_v))
{
}

/*
 * Computes the density of fdistribu
*/
void FluidMoments::operator()(double& density, DViewVx const fdistribu, FluidMoments::MomentDensity)
{
    density = m_integrate_v(fdistribu);
}

/*
 * Computes the density of allfdistribu
*/
void FluidMoments::operator()(
        DSpanSpX const density,
        DViewSpXVx const allfdistribu,
        FluidMoments::MomentDensity)
{
    ddc::for_each(
            ddc::policies::parallel_host,
            ddc::get_domain<IDimSp, IDimX>(allfdistribu),
            [&](IndexSpX const ispx) {
                (*this)(density(ispx), allfdistribu[ispx], FluidMoments::s_density);
            });
}

/*
 * Computes the mean_velocity of fdistribu, using its density
*/
void FluidMoments::operator()(
        double& mean_velocity,
        DViewVx const fdistribu,
        double const& density,
        FluidMoments::MomentVelocity)
{
    DFieldVx integrand(fdistribu.domain());
    ddc::for_each(ddc::policies::parallel_host, fdistribu.domain(), [&](IndexVx const ivx) {
        CoordVx const coordv = ddc::coordinate(ivx);
        integrand(ivx) = coordv * fdistribu(ivx);
    });

    mean_velocity = m_integrate_v(integrand) / density;
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
    ddc::for_each(
            ddc::policies::parallel_host,
            allfdistribu.domain(),
            [&](IndexSpXVx const ispxvx) {
                CoordVx const coordv = ddc::coordinate(ddc::select<IDimVx>(ispxvx));
                integrand(ispxvx) = coordv * allfdistribu(ispxvx);
            });

    ddc::for_each(
            ddc::policies::parallel_host,
            ddc::get_domain<IDimSp, IDimX>(allfdistribu),
            [&](IndexSpX const ispx) {
                mean_velocity(ispx) = m_integrate_v(integrand[ispx]) / density(ispx);
            });
}

/*
 * Computes the temperature of fdistribu, using its density and mean velocity
*/
void FluidMoments::operator()(
        double& temperature,
        DViewVx const fdistribu,
        double const& density,
        double const& mean_velocity,
        FluidMoments::MomentTemperature)
{
    DFieldVx integrand(fdistribu.domain());
    ddc::for_each(ddc::policies::parallel_host, fdistribu.domain(), [&](IndexVx const ivx) {
        double const coeff = ddc::coordinate(ddc::select<IDimVx>(ivx)) - mean_velocity;
        integrand(ivx) = coeff * coeff * fdistribu(ivx);
    });

    temperature = m_integrate_v(integrand) / density;
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
    ddc::for_each(
            ddc::policies::parallel_host,
            allfdistribu.domain(),
            [&](IndexSpXVx const ispxvx) {
                double const coeff = ddc::coordinate(ddc::select<IDimVx>(ispxvx))
                                     - mean_velocity(ddc::select<IDimSp, IDimX>(ispxvx));
                integrand(ispxvx) = coeff * coeff * allfdistribu(ispxvx);
            });

    ddc::for_each(
            ddc::policies::parallel_host,
            ddc::get_domain<IDimSp, IDimX>(allfdistribu),
            [&](IndexSpX const ispx) {
                temperature(ispx) = m_integrate_v(integrand[ispx]) / density(ispx);
            });
}
