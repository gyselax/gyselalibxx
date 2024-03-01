#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

/**
* @brief Compute the collisionality spatial profile.
* @param[inout] nustar_profile The collisionality profile.
* @param[in] nustar0 normalized collisionality coefficient.
*/
void compute_nustar_profile(DSpanSpX nustar_profile, double nustar0);

/**
* @brief Compute the collision frequency for each species.
* @param[inout] collfreq A Span representing the collision frequency for each species.
* @param[in] nustar_profile The collisionality profile.
* @param[in] density The density of each species.
* @param[in] temperature The temperature of each species.
*/
void compute_collfreq(
        DSpanSpX collfreq,
        DViewSpX nustar_profile,
        DViewSpX density,
        DViewSpX temperature);

/**
* @brief Compute the intra species collision operator diffusion coefficient.
* @param[inout] Dcoll A Span representing the diffusion coefficient.
* @param[in] collfreq The collision frequency for each species.
* @param[in] density The density of each species.
* @param[in] temperature The temperature of each species.
*/
template <class IDimension>
void compute_Dcoll(
        ddc::ChunkSpan<
                double,
                ddc::DiscreteDomain<IDimSp, IDimX, IDimension>,
                std::experimental::layout_right,
                Kokkos::DefaultExecutionSpace::memory_space> Dcoll,
        DViewSpX collfreq,
        DViewSpX density,
        DViewSpX temperature)
{
    ddc::for_each(
            ddc::policies::parallel_device,
            Dcoll.domain(),
            KOKKOS_LAMBDA(ddc::DiscreteElement<IDimSp, IDimX, IDimension> const ispxdimx) {
                double const vT(
                        Kokkos::sqrt(2. * temperature(ddc::select<IDimSp, IDimX>(ispxdimx))));
                double const v_norm(
                        Kokkos::fabs(ddc::coordinate(ddc::select<IDimension>(ispxdimx))) / vT);
                double const tol = 1.e-15;
                if (v_norm > tol) {
                    double const coeff(2. / Kokkos::sqrt(M_PI));
                    double const AD(
                            3. * Kokkos::sqrt(2. * M_PI) / 4.
                            * temperature(ddc::select<IDimSp, IDimX>(ispxdimx))
                            * collfreq(ddc::select<IDimSp, IDimX>(ispxdimx)));
                    double const inv_v_norm(1. / v_norm);
                    double const phi(Kokkos::erf(v_norm));
                    double const phi_prime(coeff * Kokkos::exp(-v_norm * v_norm));
                    double const psi((phi - v_norm * phi_prime) * 0.5 * inv_v_norm * inv_v_norm);

                    Dcoll(ispxdimx) = AD * (phi - psi) * inv_v_norm;

                } else {
                    Dcoll(ispxdimx) = Kokkos::sqrt(2)
                                      * temperature(ddc::select<IDimSp, IDimX>(ispxdimx))
                                      * collfreq(ddc::select<IDimSp, IDimX>(ispxdimx));
                }
            });
}

/**
* @brief Compute the velocity derivative of the collision operator diffusion coefficient.
* @param[inout] dvDcoll A Span representing the derivative of the diffusion coefficient.
* @param[in] collfreq The collision frequency for each species.
* @param[in] density The density of each species.
* @param[in] temperature The temperature of each species.
*/
template <class IDimension>
void compute_dvDcoll(
        ddc::ChunkSpan<
                double,
                ddc::DiscreteDomain<IDimSp, IDimX, IDimension>,
                std::experimental::layout_right,
                Kokkos::DefaultExecutionSpace::memory_space> dvDcoll,
        DViewSpX collfreq,
        DViewSpX density,
        DViewSpX temperature)
{
    ddc::for_each(
            ddc::policies::parallel_device,
            dvDcoll.domain(),
            KOKKOS_LAMBDA(ddc::DiscreteElement<IDimSp, IDimX, IDimension> const ispxdimx) {
                double const vT(
                        Kokkos::sqrt(2. * temperature(ddc::select<IDimSp, IDimX>(ispxdimx))));
                double const v_norm(
                        Kokkos::fabs(ddc::coordinate(ddc::select<IDimension>(ispxdimx))) / vT);
                double const tol = 1.e-15;
                if (v_norm > tol) {
                    double const coeff(2. / Kokkos::sqrt(M_PI));
                    double const AD(
                            3. * Kokkos::sqrt(2. * M_PI) / 4.
                            * temperature(ddc::select<IDimSp, IDimX>(ispxdimx))
                            * collfreq(ddc::select<IDimSp, IDimX>(ispxdimx)));
                    double const inv_v_norm(1. / v_norm);
                    double const phi(Kokkos::erf(v_norm));
                    double const phi_prime(coeff * Kokkos::exp(-v_norm * v_norm));
                    double const psi((phi - v_norm * phi_prime) * 0.5 * inv_v_norm * inv_v_norm);

                    double const sign(
                            ddc::coordinate(ddc::select<IDimension>(ispxdimx))
                            / Kokkos::fabs(ddc::coordinate(ddc::select<IDimension>(ispxdimx))));

                    dvDcoll(ispxdimx)
                            = sign * AD
                              / Kokkos::sqrt(2 * temperature(ddc::select<IDimSp, IDimX>(ispxdimx)))
                              * inv_v_norm * inv_v_norm * (3 * psi - phi);

                } else {
                    dvDcoll(ispxdimx) = 0.;
                }
            });
}

/**
* @brief Compute the Vcoll and Tcoll coefficients, used for building the linear system.
* 
* Computation of Vcoll and Tcoll, which are the moments
* of the kernel maxwellian function of the intra species collision operator.
* Vcoll and Tcoll are defined as follows:
*  - Tcoll = Pcoll^{-1}[Imean0*Imean2 - Imean1*Imean1]
*  - Vcoll = Pcoll^{-1}[Imean4*Imean1 - Imean3*Imean2]
*  - Pcoll = Imean0*Imean4 - Imean1*Imean3
*  where the 5 integrals are defined as:
*     Imean0=<Dcoll> ;
*     Imean1=<v*Dcoll> ;
*     Imean2=<v^2*Dcoll> ;
*     Imean3=<d/dv(Dcoll)>
*     Imean4=<d/dv(v*Dcoll)>
*  The brackets <.> represent the integral in velocity: <.> = \int . dv
* 
* @param[inout] Vcoll The Vcoll coefficient.
* @param[inout] Tcoll The Tcoll coefficient.
* @param[in] allfdistribu The distribution function.
* @param[in] Dcoll The collision operator diffusion coefficient.
* @param[in] dvDcoll The collision operator derivative of the diffusion coefficient.
*/
template <class IDimension>
void compute_Vcoll_Tcoll(
        DSpanSpX Vcoll,
        DSpanSpX Tcoll,
        DViewSpXVx allfdistribu,
        ddc::ChunkSpan<
                double,
                ddc::DiscreteDomain<IDimSp, IDimX, IDimension>,
                std::experimental::layout_right,
                Kokkos::DefaultExecutionSpace::memory_space> Dcoll,
        ddc::ChunkSpan<
                double,
                ddc::DiscreteDomain<IDimSp, IDimX, IDimension>,
                std::experimental::layout_right,
                Kokkos::DefaultExecutionSpace::memory_space> dvDcoll)
{
    host_t<DFieldVx> const quadrature_coeffs_host(
            trapezoid_quadrature_coefficients(ddc::get_domain<IDimVx>(allfdistribu)));
    auto quadrature_coeffs_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    auto quadrature_coeffs = quadrature_coeffs_alloc.span_view();

    // computation of the integrands
    DFieldSpXVx I0mean_integrand_alloc(allfdistribu.domain());
    DFieldSpXVx I1mean_integrand_alloc(allfdistribu.domain());
    DFieldSpXVx I2mean_integrand_alloc(allfdistribu.domain());
    DFieldSpXVx I3mean_integrand_alloc(allfdistribu.domain());
    DFieldSpXVx I4mean_integrand_alloc(allfdistribu.domain());
    auto I0mean_integrand = I0mean_integrand_alloc.span_view();
    auto I1mean_integrand = I1mean_integrand_alloc.span_view();
    auto I2mean_integrand = I2mean_integrand_alloc.span_view();
    auto I3mean_integrand = I3mean_integrand_alloc.span_view();
    auto I4mean_integrand = I4mean_integrand_alloc.span_view();

    ddc::for_each(
            ddc::policies::parallel_device,
            allfdistribu.domain(),
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                ddc::DiscreteElement<IDimension> const idimx(ddc::select<IDimVx>(ispxvx).uid() + 1);
                ddc::DiscreteElement<IDimSp, IDimX, IDimension>
                        ispxdimx(ddc::select<IDimSp>(ispxvx), ddc::select<IDimX>(ispxvx), idimx);
                CoordVx const coordv = ddc::coordinate(ddc::select<IDimVx>(ispxvx));
                I0mean_integrand(ispxvx) = Dcoll(ispxdimx) * allfdistribu(ispxvx);
                I1mean_integrand(ispxvx) = I0mean_integrand(ispxvx) * coordv;
                I2mean_integrand(ispxvx) = I1mean_integrand(ispxvx) * coordv;
                I3mean_integrand(ispxvx) = dvDcoll(ispxdimx) * allfdistribu(ispxvx);
                I4mean_integrand(ispxvx)
                        = I0mean_integrand(ispxvx) + I3mean_integrand(ispxvx) * coordv;
            });


    // computation of the integrals over the Vx direction
    IDomainSpX grid_sp_x(allfdistribu.domain<IDimSp, IDimX>());
    DFieldSpX I0mean_alloc(grid_sp_x);
    DFieldSpX I1mean_alloc(grid_sp_x);
    DFieldSpX I2mean_alloc(grid_sp_x);
    DFieldSpX I3mean_alloc(grid_sp_x);
    DFieldSpX I4mean_alloc(grid_sp_x);
    auto I0mean = I0mean_alloc.span_view();
    auto I1mean = I1mean_alloc.span_view();
    auto I2mean = I2mean_alloc.span_view();
    auto I3mean = I3mean_alloc.span_view();
    auto I4mean = I4mean_alloc.span_view();
    ddc::fill(I0mean, 0.);
    ddc::fill(I1mean, 0.);
    ddc::fill(I2mean, 0.);
    ddc::fill(I3mean, 0.);
    ddc::fill(I4mean, 0.);

    ddc::for_each(
            ddc::policies::parallel_device,
            grid_sp_x,
            KOKKOS_LAMBDA(IndexSpX const ispx) {
                for (IndexVx const ivx : allfdistribu.domain<IDimVx>()) {
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

/**
* @brief Compute the intra species collision operator advection coefficient.
* @param[inout] Nucoll A Span representing the advection coefficient.
* @param[in] Dcoll A Span representing the diffusion coefficient.
* @param[in] Vcoll The Vcoll coefficient.
* @param[in] Tcoll The Tcoll coefficient.
*/
template <class IDimension>
void compute_Nucoll(
        ddc::ChunkSpan<
                double,
                ddc::DiscreteDomain<IDimSp, IDimX, IDimension>,
                std::experimental::layout_right,
                Kokkos::DefaultExecutionSpace::memory_space> Nucoll,
        ddc::ChunkSpan<
                double,
                ddc::DiscreteDomain<IDimSp, IDimX, IDimension>,
                std::experimental::layout_right,
                Kokkos::DefaultExecutionSpace::memory_space> Dcoll,
        DViewSpX Vcoll,
        DViewSpX Tcoll)
{
    ddc::for_each(
            ddc::policies::parallel_device,
            Dcoll.domain(),
            KOKKOS_LAMBDA(ddc::DiscreteElement<IDimSp, IDimX, IDimension> const ispxdimx) {
                double const coordv(ddc::coordinate(ddc::select<IDimension>(ispxdimx)));
                Nucoll(ispxdimx) = -Dcoll(ispxdimx)
                                   * (coordv - Vcoll(ddc::select<IDimSp, IDimX>(ispxdimx)))
                                   / Tcoll(ddc::select<IDimSp, IDimX>(ispxdimx));
            });
}

/**
* @brief Compute the collision frequency between species a and b.
* @param[inout] collfreq_ab The collision frequency between species a and b.
* @param[in] nustar_profile The collisionality profile.
* @param[in] density The density of each species.
* @param[in] temperature The temperature of each species.
*/
void compute_collfreq_ab(
        DSpanSpX collfreq_ab,
        DViewSpX nustar_profile,
        DViewSpX density,
        DViewSpX temperature);

/**
* @brief Compute the momentum and energy exchange terms between species a and b.
* @param[inout] momentum_exchange_ab The momentum exchange term between species a and b.
* @param[inout] energy_exchange_ab The energy exchange term between species a and b.
* @param[in] collfreq_ab The collision frequency between species a and b.
* @param[in] density The density of each species.
* @param[in] mean_velocity The mean velocity of each species.
* @param[in] temperature The temperature of each species.
*/
void compute_momentum_energy_exchange(
        DSpanSpX momentum_exchange_ab,
        DSpanSpX energy_exchange_ab,
        DViewSpX collfreq_ab,
        DViewSpX density,
        DViewSpX mean_velocity,
        DViewSpX temperature);
