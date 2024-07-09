// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <fluid_moments.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

FluidMoments::FluidMoments(Quadrature<Kokkos::DefaultHostExecutionSpace, IDimVx> integrate_v)
    : m_integrate_v(integrate_v)
{
}

/*
 * Computes the density of fdistribu
*/
void FluidMoments::operator()(
        double& density,
        host_t<DViewVx> const fdistribu,
        FluidMoments::MomentDensity)
{
    density = m_integrate_v(fdistribu);
}

/*
 * Computes the density of allfdistribu
*/
void FluidMoments::operator()(
        host_t<DSpanSpX> const density,
        host_t<DViewSpXVx> const allfdistribu,
        FluidMoments::MomentDensity)
{
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        (*this)(density(ispx), allfdistribu[ispx], FluidMoments::s_density);
    });
}

/*
 * Computes the mean_velocity of fdistribu, using its density
*/
void FluidMoments::operator()(
        double& mean_velocity,
        host_t<DViewVx> const fdistribu,
        double density,
        FluidMoments::MomentVelocity)
{
    host_t<DFieldVx> integrand(fdistribu.domain());
    ddc::for_each(fdistribu.domain(), [&](IndexVx const ivx) {
        CoordVx const coordv = ddc::coordinate(ivx);
        integrand(ivx) = coordv * fdistribu(ivx);
    });

    mean_velocity = m_integrate_v(integrand) / density;
}
/*
 * Computes the mean_velocity of allfdistribu, using its density
*/
void FluidMoments::operator()(
        host_t<DSpanSpX> const mean_velocity,
        host_t<DViewSpXVx> const allfdistribu,
        host_t<DViewSpX> const density,
        FluidMoments::MomentVelocity)
{
    host_t<DFieldSpXVx> integrand(allfdistribu.domain());
    ddc::for_each(allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        CoordVx const coordv = ddc::coordinate(ddc::select<IDimVx>(ispxvx));
        integrand(ispxvx) = coordv * allfdistribu(ispxvx);
    });

    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        mean_velocity(ispx) = m_integrate_v(integrand[ispx]) / density(ispx);
    });
}

/*
 * Computes the temperature of fdistribu, using its density and mean velocity
*/
void FluidMoments::operator()(
        double& temperature,
        host_t<DViewVx> const fdistribu,
        double density,
        double mean_velocity,
        FluidMoments::MomentTemperature)
{
    host_t<DFieldVx> integrand(fdistribu.domain());
    ddc::for_each(fdistribu.domain(), [&](IndexVx const ivx) {
        double const coeff = ddc::coordinate(ddc::select<IDimVx>(ivx)) - mean_velocity;
        integrand(ivx) = coeff * coeff * fdistribu(ivx);
    });

    temperature = m_integrate_v(integrand) / density;
}
/*
 * Computes the temperature of allfdistribu, using its density and mean velocity
*/
void FluidMoments::operator()(
        host_t<DSpanSpX> const temperature,
        host_t<DViewSpXVx> const allfdistribu,
        host_t<DViewSpX> const density,
        host_t<DViewSpX> const mean_velocity,
        FluidMoments::MomentTemperature)
{
    host_t<DFieldSpXVx> integrand(allfdistribu.domain());
    ddc::for_each(allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        double const coeff = ddc::coordinate(ddc::select<IDimVx>(ispxvx))
                             - mean_velocity(ddc::select<IDimSp, IDimX>(ispxvx));
        integrand(ispxvx) = coeff * coeff * allfdistribu(ispxvx);
    });

    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        temperature(ispx) = m_integrate_v(integrand[ispx]) / density(ispx);
    });
}
