// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <fluid_moments.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

FluidMoments::FluidMoments(Quadrature<IDimVx> integrate_v) : m_integrate_v(integrate_v) {}

/*
 * Computes the density of fdistribu
*/
void FluidMoments::operator()(double& density, DViewVx const fdistribu, FluidMoments::MomentDensity)
{
    auto fdistribu_host = ddc::create_mirror_view_and_copy(fdistribu);
    density = m_integrate_v(fdistribu_host);
}

/*
 * Computes the density of allfdistribu
*/
void FluidMoments::operator()(
        DSpanSpX const density,
        DViewSpXVx const allfdistribu,
        FluidMoments::MomentDensity)
{
    auto density_host = ddc::create_mirror_view_and_copy(density);
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(allfdistribu);
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        (*this)(density_host(ispx), allfdistribu[ispx], FluidMoments::s_density);
    });
    ddc::parallel_deepcopy(density, density_host);
}

/*
 * Computes the mean_velocity of fdistribu, using its density
*/
void FluidMoments::operator()(
        double& mean_velocity,
        DViewVx const fdistribu,
        double density,
        FluidMoments::MomentVelocity)
{
    auto fdistribu_host = ddc::create_mirror_view_and_copy(fdistribu);
    host_t<DFieldVx> integrand(fdistribu.domain());

    ddc::for_each(fdistribu.domain(), [&](IndexVx const ivx) {
        CoordVx const coordv = ddc::coordinate(ivx);
        integrand(ivx) = coordv * fdistribu_host(ivx);
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
    auto mean_velocity_host = ddc::create_mirror_view_and_copy(mean_velocity);
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(allfdistribu);
    auto density_host = ddc::create_mirror_view_and_copy(density);

    host_t<DFieldSpXVx> integrand(allfdistribu.domain());
    ddc::for_each(allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        CoordVx const coordv = ddc::coordinate(ddc::select<IDimVx>(ispxvx));
        integrand(ispxvx) = coordv * allfdistribu_host(ispxvx);
    });

    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        mean_velocity_host(ispx) = m_integrate_v(integrand[ispx]) / density_host(ispx);
    });
    ddc::parallel_deepcopy(mean_velocity, mean_velocity_host);
}

/*
 * Computes the temperature of fdistribu, using its density and mean velocity
*/
void FluidMoments::operator()(
        double& temperature,
        DViewVx const fdistribu,
        double density,
        double mean_velocity,
        FluidMoments::MomentTemperature)
{
    auto fdistribu_host = ddc::create_mirror_view_and_copy(fdistribu);
    host_t<DFieldVx> integrand(fdistribu.domain());
    ddc::for_each(fdistribu.domain(), [&](IndexVx const ivx) {
        double const coeff = ddc::coordinate(ddc::select<IDimVx>(ivx)) - mean_velocity;
        integrand(ivx) = coeff * coeff * fdistribu_host(ivx);
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
    auto temperature_host = ddc::create_mirror_view_and_copy(temperature);
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(allfdistribu);
    auto density_host = ddc::create_mirror_view_and_copy(density);
    auto mean_velocity_host = ddc::create_mirror_view_and_copy(mean_velocity);

    host_t<DFieldSpXVx> integrand(allfdistribu.domain());
    ddc::for_each(allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        double const coeff = ddc::coordinate(ddc::select<IDimVx>(ispxvx))
                             - mean_velocity_host(ddc::select<IDimSp, IDimX>(ispxvx));
        integrand(ispxvx) = coeff * coeff * allfdistribu_host(ispxvx);
    });

    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        temperature_host(ispx) = m_integrate_v(integrand[ispx]) / density_host(ispx);
    });
    ddc::parallel_deepcopy(temperature, temperature_host);
}
