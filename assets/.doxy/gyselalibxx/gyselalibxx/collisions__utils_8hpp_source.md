

# File collisions\_utils.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**collisions\_utils.hpp**](collisions__utils_8hpp.md)

[Go to the documentation of this file](collisions__utils_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

void compute_nustar_profile(DFieldSpX nustar_profile, double nustar0);

void compute_collfreq(
        DFieldSpX collfreq,
        DConstFieldSpX nustar_profile,
        DConstFieldSpX density,
        DConstFieldSpX temperature);

template <class LocalGridVx>
void compute_Dcoll(
        DField<IdxRange<Species, GridX, LocalGridVx>> Dcoll,
        DConstFieldSpX collfreq,
        DConstFieldSpX density,
        DConstFieldSpX temperature)
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(Dcoll),
            KOKKOS_LAMBDA(Idx<Species, GridX, LocalGridVx> const ispxdimvx) {
                double const vT(
                        Kokkos::sqrt(2. * temperature(ddc::select<Species, GridX>(ispxdimvx))));
                double const v_norm(
                        Kokkos::fabs(ddc::coordinate(ddc::select<LocalGridVx>(ispxdimvx))) / vT);
                double const tol = 1.e-15;
                if (v_norm > tol) {
                    double const coeff(2. / Kokkos::sqrt(M_PI));
                    double const AD(
                            3. * Kokkos::sqrt(2. * M_PI) / 4.
                            * temperature(ddc::select<Species, GridX>(ispxdimvx))
                            * collfreq(ddc::select<Species, GridX>(ispxdimvx)));
                    double const inv_v_norm(1. / v_norm);
                    double const phi(Kokkos::erf(v_norm));
                    double const phi_prime(coeff * Kokkos::exp(-v_norm * v_norm));
                    double const psi((phi - v_norm * phi_prime) * 0.5 * inv_v_norm * inv_v_norm);

                    Dcoll(ispxdimvx) = AD * (phi - psi) * inv_v_norm;

                } else {
                    Dcoll(ispxdimvx) = Kokkos::sqrt(2)
                                       * temperature(ddc::select<Species, GridX>(ispxdimvx))
                                       * collfreq(ddc::select<Species, GridX>(ispxdimvx));
                }
            });
}

template <class LocalGridVx>
void compute_dvDcoll(
        DField<IdxRange<Species, GridX, LocalGridVx>> dvDcoll,
        DConstFieldSpX collfreq,
        DConstFieldSpX density,
        DConstFieldSpX temperature)
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(dvDcoll),
            KOKKOS_LAMBDA(Idx<Species, GridX, LocalGridVx> const ispxdimvx) {
                double const vT(
                        Kokkos::sqrt(2. * temperature(ddc::select<Species, GridX>(ispxdimvx))));
                double const v_norm(
                        Kokkos::fabs(ddc::coordinate(ddc::select<LocalGridVx>(ispxdimvx))) / vT);
                double const tol = 1.e-15;
                if (v_norm > tol) {
                    double const coeff(2. / Kokkos::sqrt(M_PI));
                    double const AD(
                            3. * Kokkos::sqrt(2. * M_PI) / 4.
                            * temperature(ddc::select<Species, GridX>(ispxdimvx))
                            * collfreq(ddc::select<Species, GridX>(ispxdimvx)));
                    double const inv_v_norm(1. / v_norm);
                    double const phi(Kokkos::erf(v_norm));
                    double const phi_prime(coeff * Kokkos::exp(-v_norm * v_norm));
                    double const psi((phi - v_norm * phi_prime) * 0.5 * inv_v_norm * inv_v_norm);

                    double const sign(
                            ddc::coordinate(ddc::select<LocalGridVx>(ispxdimvx))
                            / Kokkos::fabs(ddc::coordinate(ddc::select<LocalGridVx>(ispxdimvx))));

                    dvDcoll(ispxdimvx)
                            = sign * AD
                              / Kokkos::sqrt(
                                      2 * temperature(ddc::select<Species, GridX>(ispxdimvx)))
                              * inv_v_norm * inv_v_norm * (3 * psi - phi);

                } else {
                    dvDcoll(ispxdimvx) = 0.;
                }
            });
}

template <class LocalGridVx>
void compute_Vcoll_Tcoll(
        DFieldSpX Vcoll,
        DFieldSpX Tcoll,
        DConstFieldSpXVx allfdistribu,
        DField<IdxRange<Species, GridX, LocalGridVx>> Dcoll,
        DField<IdxRange<Species, GridX, LocalGridVx>> dvDcoll)
{
    DFieldMemVx const quadrature_coeffs_alloc(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(
                    get_idx_range<GridVx>(allfdistribu)));
    DConstFieldVx const quadrature_coeffs = get_const_field(quadrature_coeffs_alloc);

    // computation of the integrands
    DFieldMemSpXVx I0mean_integrand_alloc(get_idx_range(allfdistribu));
    DFieldMemSpXVx I1mean_integrand_alloc(get_idx_range(allfdistribu));
    DFieldMemSpXVx I2mean_integrand_alloc(get_idx_range(allfdistribu));
    DFieldMemSpXVx I3mean_integrand_alloc(get_idx_range(allfdistribu));
    DFieldMemSpXVx I4mean_integrand_alloc(get_idx_range(allfdistribu));
    auto I0mean_integrand = get_field(I0mean_integrand_alloc);
    auto I1mean_integrand = get_field(I1mean_integrand_alloc);
    auto I2mean_integrand = get_field(I2mean_integrand_alloc);
    auto I3mean_integrand = get_field(I3mean_integrand_alloc);
    auto I4mean_integrand = get_field(I4mean_integrand_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(allfdistribu),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                Idx<LocalGridVx> const idimvx(ddc::select<GridVx>(ispxvx).uid() + 1);
                Idx<Species, GridX, LocalGridVx>
                        ispxdimvx(ddc::select<Species>(ispxvx), ddc::select<GridX>(ispxvx), idimvx);
                CoordVx const coordv = ddc::coordinate(ddc::select<GridVx>(ispxvx));
                I0mean_integrand(ispxvx) = Dcoll(ispxdimvx) * allfdistribu(ispxvx);
                I1mean_integrand(ispxvx) = I0mean_integrand(ispxvx) * coordv;
                I2mean_integrand(ispxvx) = I1mean_integrand(ispxvx) * coordv;
                I3mean_integrand(ispxvx) = dvDcoll(ispxdimvx) * allfdistribu(ispxvx);
                I4mean_integrand(ispxvx)
                        = I0mean_integrand(ispxvx) + I3mean_integrand(ispxvx) * coordv;
            });


    // computation of the integrals over the Vx direction
    IdxRangeSpX grid_sp_x(get_idx_range<Species, GridX>(allfdistribu));
    DFieldMemSpX I0mean_alloc(grid_sp_x);
    DFieldMemSpX I1mean_alloc(grid_sp_x);
    DFieldMemSpX I2mean_alloc(grid_sp_x);
    DFieldMemSpX I3mean_alloc(grid_sp_x);
    DFieldMemSpX I4mean_alloc(grid_sp_x);
    auto I0mean = get_field(I0mean_alloc);
    auto I1mean = get_field(I1mean_alloc);
    auto I2mean = get_field(I2mean_alloc);
    auto I3mean = get_field(I3mean_alloc);
    auto I4mean = get_field(I4mean_alloc);
    ddc::parallel_fill(I0mean, 0.);
    ddc::parallel_fill(I1mean, 0.);
    ddc::parallel_fill(I2mean, 0.);
    ddc::parallel_fill(I3mean, 0.);
    ddc::parallel_fill(I4mean, 0.);

    IdxRangeVx const idx_range_vx(get_idx_range<GridVx>(allfdistribu));

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            grid_sp_x,
            KOKKOS_LAMBDA(IdxSpX const ispx) {
                for (IdxVx const ivx : idx_range_vx) {
                    I0mean(ispx) += quadrature_coeffs(ivx) * I0mean_integrand(ispx, ivx);
                    I1mean(ispx) += quadrature_coeffs(ivx) * I1mean_integrand(ispx, ivx);
                    I2mean(ispx) += quadrature_coeffs(ivx) * I2mean_integrand(ispx, ivx);
                    I3mean(ispx) += quadrature_coeffs(ivx) * I3mean_integrand(ispx, ivx);
                    I4mean(ispx) += quadrature_coeffs(ivx) * I4mean_integrand(ispx, ivx);
                }

                double const inv_Pcoll(
                        1. / (I0mean(ispx) * I4mean(ispx) - I1mean(ispx) * I3mean(ispx)));
                Vcoll(ispx)
                        = inv_Pcoll * (I1mean(ispx) * I4mean(ispx) - I2mean(ispx) * I3mean(ispx));
                Tcoll(ispx)
                        = inv_Pcoll * (I0mean(ispx) * I2mean(ispx) - I1mean(ispx) * I1mean(ispx));
            });
}

template <class LocalGridVx>
void compute_Nucoll(
        DField<IdxRange<Species, GridX, LocalGridVx>> Nucoll,
        DField<IdxRange<Species, GridX, LocalGridVx>> Dcoll,
        DConstFieldSpX Vcoll,
        DConstFieldSpX Tcoll)
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(Dcoll),
            KOKKOS_LAMBDA(Idx<Species, GridX, LocalGridVx> const ispxdimvx) {
                double const coordv(ddc::coordinate(ddc::select<LocalGridVx>(ispxdimvx)));
                Nucoll(ispxdimvx) = -Dcoll(ispxdimvx)
                                    * (coordv - Vcoll(ddc::select<Species, GridX>(ispxdimvx)))
                                    / Tcoll(ddc::select<Species, GridX>(ispxdimvx));
            });
}

void compute_collfreq_ab(
        DFieldSpX collfreq_ab,
        DConstFieldSpX nustar_profile,
        DConstFieldSpX density,
        DConstFieldSpX temperature);

void compute_momentum_energy_exchange(
        DFieldSpX momentum_exchange_ab,
        DFieldSpX energy_exchange_ab,
        DConstFieldSpX collfreq_ab,
        DConstFieldSpX density,
        DConstFieldSpX mean_velocity,
        DConstFieldSpX temperature);
```


