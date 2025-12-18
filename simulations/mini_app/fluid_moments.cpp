// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "fluid_moments.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

FluidMoments::FluidMoments(
        Quadrature<IdxRangeVparMu, IdxRangeSpTor3DV2D, Kokkos::DefaultExecutionSpace::memory_space>
                integrate_vparmu)
    : m_integrate_vparmu(integrate_vparmu)
{
}

/*
 * Computes the density of allfdistribu by integrating over vpar and mu
*/
void FluidMoments::operator()(
        DField<IdxRangeSpTor3D> const density,
        DConstField<IdxRangeSpTor3DV2D> const allfdistribu,
        FluidMoments::MomentDensity)
{
    m_integrate_vparmu(Kokkos::DefaultExecutionSpace(), density, allfdistribu);
}

/*
 * Computes the mean_velocity (vpar) of allfdistribu, using its density
*/
void FluidMoments::operator()(
        DField<IdxRangeSpTor3D> const mean_velocity,
        DConstField<IdxRangeSpTor3DV2D> const allfdistribu,
        DConstField<IdxRangeSpTor3D> const density,
        FluidMoments::MomentVelocity)
{
    m_integrate_vparmu(
            Kokkos::DefaultExecutionSpace(),
            mean_velocity,
            KOKKOS_LAMBDA(IdxSpTor3DV2D const ispgrid) {
                IdxSpTor3D isptor123(ispgrid);
                IdxVpar ivpar(ddc::select<GridVpar>(ispgrid));
                CoordVpar const coordvpar = ddc::coordinate(ivpar);
                return coordvpar * allfdistribu(ispgrid) / density(isptor123);
            });
}

/*
 * Computes the temperature of allfdistribu, using its density and mean velocity
 * Temperature is computed as: T = <(vpar - Upar)^2> where Upar is the mean velocity
*/
void FluidMoments::operator()(
        DField<IdxRangeSpTor3D> const temperature,
        DConstField<IdxRangeSpTor3DV2D> const allfdistribu,
        DConstField<IdxRangeSpTor3D> const density,
        DConstField<IdxRangeSpTor3D> const mean_velocity,
        FluidMoments::MomentTemperature)
{
    m_integrate_vparmu(
            Kokkos::DefaultExecutionSpace(),
            temperature,
            KOKKOS_LAMBDA(IdxSpTor3DV2D const ispgrid) {
                IdxSpTor3D isptor123(ispgrid);
                IdxVpar ivpar(ddc::select<GridVpar>(ispgrid));
                double const coeff = ddc::coordinate(ivpar) - mean_velocity(isptor123);
                return coeff * coeff * allfdistribu(ispgrid) / density(isptor123);
            });
}
