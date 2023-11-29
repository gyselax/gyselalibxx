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
        ddc::ChunkSpan<double, ddc::DiscreteDomain<IDimSp, IDimX, IDimension>> Dcoll,
        DViewSpX collfreq,
        DViewSpX density,
        DViewSpX temperature)
{
    ddc::for_each(Dcoll.domain(), [&](auto const ispxvx) {
        double const vT(std::sqrt(2. * temperature(ddc::select<IDimSp, IDimX>(ispxvx))));
        double const v_norm(std::fabs(ddc::coordinate(ddc::select<IDimension>(ispxvx))) / vT);
        double const tol = 1.e-15;
        if (v_norm > tol) {
            double const coeff(2. / std::sqrt(M_PI));
            double const AD(
                    3. * std::sqrt(2. * M_PI) / 4. * temperature(ddc::select<IDimSp, IDimX>(ispxvx))
                    * collfreq(ddc::select<IDimSp, IDimX>(ispxvx)));
            double const inv_v_norm(1. / v_norm);
            double const phi(std::erf(v_norm));
            double const phi_prime(coeff * std::exp(-v_norm * v_norm));
            double const psi((phi - v_norm * phi_prime) * 0.5 * inv_v_norm * inv_v_norm);

            Dcoll(ispxvx) = AD * (phi - psi) * inv_v_norm;

        } else {
            Dcoll(ispxvx) = std::sqrt(2) * temperature(ddc::select<IDimSp, IDimX>(ispxvx))
                            * collfreq(ddc::select<IDimSp, IDimX>(ispxvx));
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
        ddc::ChunkSpan<double, ddc::DiscreteDomain<IDimSp, IDimX, IDimension>> dvDcoll,
        DViewSpX collfreq,
        DViewSpX density,
        DViewSpX temperature)
{
    ddc::for_each(dvDcoll.domain(), [&](auto const ispxvx) {
        double const vT(std::sqrt(2. * temperature(ddc::select<IDimSp, IDimX>(ispxvx))));
        double const v_norm(std::fabs(ddc::coordinate(ddc::select<IDimension>(ispxvx))) / vT);
        double const tol = 1.e-15;
        if (v_norm > tol) {
            double const coeff(2. / std::sqrt(M_PI));
            double const AD(
                    3. * std::sqrt(2. * M_PI) / 4. * temperature(ddc::select<IDimSp, IDimX>(ispxvx))
                    * collfreq(ddc::select<IDimSp, IDimX>(ispxvx)));
            double const inv_v_norm(1. / v_norm);
            double const phi(std::erf(v_norm));
            double const phi_prime(coeff * std::exp(-v_norm * v_norm));
            double const psi((phi - v_norm * phi_prime) * 0.5 * inv_v_norm * inv_v_norm);

            double const sign(
                    ddc::coordinate(ddc::select<IDimension>(ispxvx))
                    / std::fabs(ddc::coordinate(ddc::select<IDimension>(ispxvx))));

            dvDcoll(ispxvx) = sign * AD
                              / std::sqrt(2 * temperature(ddc::select<IDimSp, IDimX>(ispxvx)))
                              * inv_v_norm * inv_v_norm * (3 * psi - phi);

        } else {
            dvDcoll(ispxvx) = 0.;
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
        ddc::ChunkSpan<double const, ddc::DiscreteDomain<IDimSp, IDimX, IDimension>> Dcoll,
        ddc::ChunkSpan<double const, ddc::DiscreteDomain<IDimSp, IDimX, IDimension>> dvDcoll)
{
    DFieldVx const quadrature_coeffs
            = trapezoid_quadrature_coefficients(ddc::get_domain<IDimVx>(allfdistribu));
    Quadrature<IDimVx> const integrate_v(quadrature_coeffs);

    // computation of the integrands
    DFieldSpXVx I0mean_integrand(allfdistribu.domain());
    DFieldSpXVx I1mean_integrand(allfdistribu.domain());
    DFieldSpXVx I2mean_integrand(allfdistribu.domain());
    DFieldSpXVx I3mean_integrand(allfdistribu.domain());
    DFieldSpXVx I4mean_integrand(allfdistribu.domain());
    ddc::for_each(allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        ddc::DiscreteElement<IDimension> const idimx(ddc::select<IDimVx>(ispxvx).uid() + 1);
        ddc::DiscreteElement<IDimSp, IDimX, IDimension>
                ispxdimx(ddc::select<IDimSp>(ispxvx), ddc::select<IDimX>(ispxvx), idimx);
        CoordVx const coordv = ddc::coordinate(ddc::select<IDimVx>(ispxvx));
        I0mean_integrand(ispxvx) = Dcoll(ispxdimx) * allfdistribu(ispxvx);
        I1mean_integrand(ispxvx) = I0mean_integrand(ispxvx) * coordv;
        I2mean_integrand(ispxvx) = I1mean_integrand(ispxvx) * coordv;
        I3mean_integrand(ispxvx) = dvDcoll(ispxdimx) * allfdistribu(ispxvx);
        I4mean_integrand(ispxvx) = I0mean_integrand(ispxvx) + I3mean_integrand(ispxvx) * coordv;
    });

    // computation of the integrals over the Vx direction
    DFieldSpX I0mean(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX I1mean(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX I2mean(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX I3mean(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX I4mean(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        I0mean(ispx) = integrate_v(I0mean_integrand[ispx]);
        I1mean(ispx) = integrate_v(I1mean_integrand[ispx]);
        I2mean(ispx) = integrate_v(I2mean_integrand[ispx]);
        I3mean(ispx) = integrate_v(I3mean_integrand[ispx]);
        I4mean(ispx) = integrate_v(I4mean_integrand[ispx]);

        double const inv_Pcoll(1. / (I0mean(ispx) * I4mean(ispx) - I1mean(ispx) * I3mean(ispx)));
        Vcoll(ispx) = inv_Pcoll * (I1mean(ispx) * I4mean(ispx) - I2mean(ispx) * I3mean(ispx));
        Tcoll(ispx) = inv_Pcoll * (I0mean(ispx) * I2mean(ispx) - I1mean(ispx) * I1mean(ispx));
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
        ddc::ChunkSpan<double, ddc::DiscreteDomain<IDimSp, IDimX, IDimension>> Nucoll,
        ddc::ChunkSpan<double const, ddc::DiscreteDomain<IDimSp, IDimX, IDimension>> Dcoll,
        DViewSpX Vcoll,
        DViewSpX Tcoll)
{
    ddc::for_each(Dcoll.domain(), [&](auto const ispxdimx) {
        double const coordv(ddc::coordinate(ddc::select<IDimension>(ispxdimx)));
        Nucoll(ispxdimx) = -Dcoll(ispxdimx) * (coordv - Vcoll(ddc::select<IDimSp, IDimX>(ispxdimx)))
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
        DSpanSp collfreq_ab,
        DViewSp nustar_profile,
        DViewSp density,
        DViewSp temperature);

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
        DSpanSp momentum_exchange_ab,
        DSpanSp energy_exchange_ab,
        DViewSp collfreq_ab,
        DViewSp density,
        DViewSp mean_velocity,
        DViewSp temperature);
