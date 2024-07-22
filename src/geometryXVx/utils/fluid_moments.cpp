// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <fluid_moments.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

FluidMoments::FluidMoments(
        Quadrature<IDomainVx, IDomainSpXVx, Kokkos::DefaultExecutionSpace::memory_space>
                integrate_v)
    : m_integrate_v(integrate_v)
{
}

/*
 * Computes the density of fdistribu
*/
void FluidMoments::operator()(double& density, DViewVx const fdistribu, FluidMoments::MomentDensity)
{
    density = m_integrate_v(Kokkos::DefaultExecutionSpace(), fdistribu);
}

/*
 * Computes the density of allfdistribu
*/
void FluidMoments::operator()(
        DSpanSpX const density,
        DViewSpXVx const allfdistribu,
        FluidMoments::MomentDensity)
{
    m_integrate_v(Kokkos::DefaultExecutionSpace(), density, allfdistribu);
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
    DFieldVx integrand_alloc(fdistribu.domain());
    DSpanVx integrand = integrand_alloc.span_view();

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            fdistribu.domain(),
            KOKKOS_LAMBDA(IndexVx const ivx) {
                CoordVx const coordv = ddc::coordinate(ivx);
                integrand(ivx) = coordv * fdistribu(ivx);
            });

    mean_velocity = m_integrate_v(Kokkos::DefaultExecutionSpace(), integrand) / density;
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
    m_integrate_v(
            Kokkos::DefaultExecutionSpace(),
            mean_velocity,
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                IndexSpX ispx(ispxvx);
                IndexVx ivx(ispxvx);
                CoordVx const coordv = ddc::coordinate(ddc::select<IDimVx>(ispxvx));
                return coordv * allfdistribu(ispxvx) / density(ispx);
            });
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
    temperature = m_integrate_v(
                          Kokkos::DefaultExecutionSpace(),
                          KOKKOS_LAMBDA(IndexVx const ivx) {
                              double const coeff
                                      = ddc::coordinate(ddc::select<IDimVx>(ivx)) - mean_velocity;
                              return coeff * coeff * fdistribu(ivx);
                          })
                  / density;
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
    m_integrate_v(
            Kokkos::DefaultExecutionSpace(),
            temperature,
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                IndexSpX ispx(ispxvx);
                IndexVx ivx(ispxvx);
                double const coeff = ddc::coordinate(ivx) - mean_velocity(ispx);
                return coeff * coeff * allfdistribu(ispxvx) / density(ispx);
            });
}
