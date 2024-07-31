// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <fluid_moments.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

FluidMoments::FluidMoments(
        Quadrature<IdxRangeVx, IdxRangeSpXVx, Kokkos::DefaultExecutionSpace::memory_space>
                integrate_v)
    : m_integrate_v(integrate_v)
{
}

/*
 * Computes the density of fdistribu
*/
void FluidMoments::operator()(
        double& density,
        DConstFieldVx const fdistribu,
        FluidMoments::MomentDensity)
{
    density = m_integrate_v(Kokkos::DefaultExecutionSpace(), fdistribu);
}

/*
 * Computes the density of allfdistribu
*/
void FluidMoments::operator()(
        DFieldSpX const density,
        DConstFieldSpXVx const allfdistribu,
        FluidMoments::MomentDensity)
{
    m_integrate_v(Kokkos::DefaultExecutionSpace(), density, allfdistribu);
}

/*
 * Computes the mean_velocity of fdistribu, using its density
*/
void FluidMoments::operator()(
        double& mean_velocity,
        DConstFieldVx const fdistribu,
        double density,
        FluidMoments::MomentVelocity)
{
    DFieldMemVx integrand_alloc(get_idx_range(fdistribu));
    DFieldVx integrand = get_field(integrand_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(fdistribu),
            KOKKOS_LAMBDA(IdxVx const ivx) {
                CoordVx const coordv = ddc::coordinate(ivx);
                integrand(ivx) = coordv * fdistribu(ivx);
            });

    mean_velocity = m_integrate_v(Kokkos::DefaultExecutionSpace(), integrand) / density;
}
/*
 * Computes the mean_velocity of allfdistribu, using its density
*/
void FluidMoments::operator()(
        DFieldSpX const mean_velocity,
        DConstFieldSpXVx const allfdistribu,
        DConstFieldSpX const density,
        FluidMoments::MomentVelocity)
{
    m_integrate_v(
            Kokkos::DefaultExecutionSpace(),
            mean_velocity,
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                IdxSpX ispx(ispxvx);
                IdxVx ivx(ispxvx);
                CoordVx const coordv = ddc::coordinate(ddc::select<GridVx>(ispxvx));
                return coordv * allfdistribu(ispxvx) / density(ispx);
            });
}

/*
 * Computes the temperature of fdistribu, using its density and mean velocity
*/
void FluidMoments::operator()(
        double& temperature,
        DConstFieldVx const fdistribu,
        double density,
        double mean_velocity,
        FluidMoments::MomentTemperature)
{
    temperature = m_integrate_v(
                          Kokkos::DefaultExecutionSpace(),
                          KOKKOS_LAMBDA(IdxVx const ivx) {
                              double const coeff
                                      = ddc::coordinate(ddc::select<GridVx>(ivx)) - mean_velocity;
                              return coeff * coeff * fdistribu(ivx);
                          })
                  / density;
}
/*
 * Computes the temperature of allfdistribu, using its density and mean velocity
*/
void FluidMoments::operator()(
        DFieldSpX const temperature,
        DConstFieldSpXVx const allfdistribu,
        DConstFieldSpX const density,
        DConstFieldSpX const mean_velocity,
        FluidMoments::MomentTemperature)
{
    m_integrate_v(
            Kokkos::DefaultExecutionSpace(),
            temperature,
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                IdxSpX ispx(ispxvx);
                IdxVx ivx(ispxvx);
                double const coeff = ddc::coordinate(ivx) - mean_velocity(ispx);
                return coeff * coeff * allfdistribu(ispxvx) / density(ispx);
            });
}
