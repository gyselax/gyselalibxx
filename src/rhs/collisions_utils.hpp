#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

IndexSp iion();

void compute_nustar_profile(DSpanSpX nustar_profile, double nustar0);

// intra species collision operator helper functions
void compute_collfreq(
        DSpanSpX collfreq,
        DViewSpX nustar_profile,
        DViewSpX density,
        DViewSpX temperature);

/**
 * Computes the diffusion coefficent Dcoll
 */
template <class IDimension>
void compute_Dcoll(
        ddc::ChunkSpan<double, DiscreteDomain<IDimSp, IDimX, IDimension>> Dcoll,
        DViewSpX collfreq,
        DViewSpX density,
        DViewSpX temperature)
{
    for_each(policies::parallel_host, Dcoll.domain(), [&](auto const ispxvx) {
        double const vT(std::sqrt(2. * temperature(select<IDimSp, IDimX>(ispxvx))));
        double const v_norm(std::fabs(coordinate(select<IDimension>(ispxvx))) / vT);
        double const tol = 1.e-15;
        if (v_norm > tol) {
            double const coeff(2. / std::sqrt(M_PI));
            double const AD(
                    3. * std::sqrt(2. * M_PI) / 4. * temperature(select<IDimSp, IDimX>(ispxvx))
                    * collfreq(select<IDimSp, IDimX>(ispxvx)));
            double const inv_v_norm(1. / v_norm);
            double const phi(std::erf(v_norm));
            double const phi_prime(coeff * std::exp(-v_norm * v_norm));
            double const psi((phi - v_norm * phi_prime) * 0.5 * inv_v_norm * inv_v_norm);

            Dcoll(ispxvx) = AD * (phi - psi) * inv_v_norm;

        } else {
            Dcoll(ispxvx) = std::sqrt(2) * temperature(select<IDimSp, IDimX>(ispxvx))
                            * collfreq(select<IDimSp, IDimX>(ispxvx));
        }
    });
}

template <class IDimension>
void compute_dvDcoll(
        ddc::ChunkSpan<double, DiscreteDomain<IDimSp, IDimX, IDimension>> dvDcoll,
        DViewSpX collfreq,
        DViewSpX density,
        DViewSpX temperature)
{
    for_each(policies::parallel_host, dvDcoll.domain(), [&](auto const ispxvx) {
        double const vT(std::sqrt(2. * temperature(select<IDimSp, IDimX>(ispxvx))));
        double const v_norm(std::fabs(coordinate(select<IDimension>(ispxvx))) / vT);
        double const tol = 1.e-15;
        if (v_norm > tol) {
            double const coeff(2. / std::sqrt(M_PI));
            double const AD(
                    3. * std::sqrt(2. * M_PI) / 4. * temperature(select<IDimSp, IDimX>(ispxvx))
                    * collfreq(select<IDimSp, IDimX>(ispxvx)));
            double const inv_v_norm(1. / v_norm);
            double const phi(std::erf(v_norm));
            double const phi_prime(coeff * std::exp(-v_norm * v_norm));
            double const psi((phi - v_norm * phi_prime) * 0.5 * inv_v_norm * inv_v_norm);

            double const sign(
                    coordinate(select<IDimension>(ispxvx))
                    / std::fabs(coordinate(select<IDimension>(ispxvx))));

            dvDcoll(ispxvx) = sign * AD / std::sqrt(2 * temperature(select<IDimSp, IDimX>(ispxvx)))
                              * inv_v_norm * inv_v_norm * (3 * psi - phi);

        } else {
            dvDcoll(ispxvx) = 0.;
        }
    });
}

/**
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
 */
template <class IDimension>
void compute_Vcoll_Tcoll(
        DSpanSpX Vcoll,
        DSpanSpX Tcoll,
        DViewSpXVx allfdistribu,
        ddc::ChunkSpan<double const, DiscreteDomain<IDimSp, IDimX, IDimension>> Dcoll,
        ddc::ChunkSpan<double const, DiscreteDomain<IDimSp, IDimX, IDimension>> dvDcoll)
{
    Quadrature<IDimVx> const integrate_v(
            trapezoid_quadrature_coefficients(get_domain<IDimVx>(allfdistribu)));

    // computation of the integrands
    DFieldSpXVx I0mean_integrand(allfdistribu.domain());
    DFieldSpXVx I1mean_integrand(allfdistribu.domain());
    DFieldSpXVx I2mean_integrand(allfdistribu.domain());
    DFieldSpXVx I3mean_integrand(allfdistribu.domain());
    DFieldSpXVx I4mean_integrand(allfdistribu.domain());
    for_each(policies::parallel_host, allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        DiscreteElement<IDimension> const idimx(select<IDimVx>(ispxvx).uid() + 1);
        DiscreteElement<IDimSp, IDimX, IDimension>
                ispxdimx(select<IDimSp>(ispxvx), select<IDimX>(ispxvx), idimx);
        CoordVx const coordv = coordinate(select<IDimVx>(ispxvx));
        I0mean_integrand(ispxvx) = Dcoll(ispxdimx) * allfdistribu(ispxvx);
        I1mean_integrand(ispxvx) = I0mean_integrand(ispxvx) * coordv;
        I2mean_integrand(ispxvx) = I1mean_integrand(ispxvx) * coordv;
        I3mean_integrand(ispxvx) = dvDcoll(ispxdimx) * allfdistribu(ispxvx);
        I4mean_integrand(ispxvx) = I0mean_integrand(ispxvx) + I3mean_integrand(ispxvx) * coordv;
    });

    // computation of the integrals over the Vx direction
    DFieldSpX I0mean(get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX I1mean(get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX I2mean(get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX I3mean(get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX I4mean(get_domain<IDimSp, IDimX>(allfdistribu));
    for_each(
            policies::parallel_host,
            get_domain<IDimSp, IDimX>(allfdistribu),
            [&](IndexSpX const ispx) {
                I0mean(ispx) = integrate_v(I0mean_integrand[ispx]);
                I1mean(ispx) = integrate_v(I1mean_integrand[ispx]);
                I2mean(ispx) = integrate_v(I2mean_integrand[ispx]);
                I3mean(ispx) = integrate_v(I3mean_integrand[ispx]);
                I4mean(ispx) = integrate_v(I4mean_integrand[ispx]);

                double const inv_Pcoll(
                        1. / (I0mean(ispx) * I4mean(ispx) - I1mean(ispx) * I3mean(ispx)));
                Vcoll(ispx)
                        = inv_Pcoll * (I1mean(ispx) * I4mean(ispx) - I2mean(ispx) * I3mean(ispx));
                Tcoll(ispx)
                        = inv_Pcoll * (I0mean(ispx) * I2mean(ispx) - I1mean(ispx) * I1mean(ispx));
            });
}

/**
 * Computes the convection coefficent Nucoll
 */
template <class IDimension>
void compute_Nucoll(
        ddc::ChunkSpan<double, DiscreteDomain<IDimSp, IDimX, IDimension>> Nucoll,
        ddc::ChunkSpan<double const, DiscreteDomain<IDimSp, IDimX, IDimension>> Dcoll,
        DViewSpX Vcoll,
        DViewSpX Tcoll)
{
    for_each(policies::parallel_host, Dcoll.domain(), [&](auto const ispxdimx) {
        double const coordv(coordinate(select<IDimension>(ispxdimx)));
        Nucoll(ispxdimx) = -Dcoll(ispxdimx) * (coordv - Vcoll(select<IDimSp, IDimX>(ispxdimx)))
                           / Tcoll(select<IDimSp, IDimX>(ispxdimx));
    });
}

// inter species collision operator helper functions
void compute_collfreq_ei(
        DSpanSpX collfreq_ei,
        DViewSpX nustar_profile,
        DViewSpX density,
        DViewSpX temperature);

void compute_momentum_energy_exchange(
        DSpanSpX momentum_exchange_ei,
        DSpanSpX momentum_exchange_ie,
        DSpanSpX energy_exchange_ei,
        DSpanSpX energy_exchange_ie,
        DViewX collfreq_ei,
        DViewSpX density,
        DViewSpX mean_velocity,
        DViewSpX temperature);
