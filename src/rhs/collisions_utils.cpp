
#include <cassert>
#include <cmath>

#include <ddc_helper.hpp>
#include <quadrature.hpp>
#include <species_info.hpp>
#include <trapezoid_quadrature.hpp>

#include "collisions_utils.hpp"

/**
 * Useful function for computing collision related quantities: 
 *  - kernel maxwellian moments
 */

IndexSp iion()
{
    return ielec() + 1;
}

/**
 * Computes the spatial profile of nustar, which is constant here,
 * but could be space dependent (for instance to have no collisions in some specific
 * parts of the simulation box).
 */
void compute_nustar_profile(DSpanSpX nustar_profile, double nustar0)
{
    double const Lx = ddcHelper::total_interval_length(get_domain<IDimX>(nustar_profile));
    for_each(policies::parallel_host, nustar_profile.domain(), [&](IndexSpX const ispx) {
        double const coeff = std::sqrt(mass(ielec()) / mass(select<IDimSp>(ispx)))
                             * std::pow(charge(select<IDimSp>(ispx)), 4) / Lx;
        nustar_profile(ispx) = coeff * nustar0;
    });
}

/**
 * Computes the space and species dependent collision frequency collfreq.
 */
void compute_collfreq(
        DSpanSpX collfreq,
        DViewSpX nustar_profile,
        DViewSpX density,
        DViewSpX temperature)
{
    for_each(policies::parallel_host, collfreq.domain(), [&](IndexSpX const ispx) {
        collfreq(ispx) = nustar_profile(ispx) * density(ispx) / std::pow(temperature(ispx), 1.5);
    });
}

/**
 * Computes the two species collision frequency collfreq_ei
 */
void compute_collfreq_ei(
        DSpanX collfreq_ei,
        DViewSpX nustar_profile,
        DViewSpX density,
        DViewSpX temperature)
{
    double const coeff(
            charge(iion()) * charge(iion())
            / (4 * std::sqrt(2) * charge(ielec()) * charge(ielec())));
    double const mass_ratio(mass(ielec()) / mass(iion()));
    for_each(policies::parallel_host, collfreq_ei.domain(), [&](IndexX const ix) {
        double const collfreq_elec(
                nustar_profile(ielec(), ix) * density(ielec(), ix)
                / std::pow(temperature(ielec(), ix), 1.5));

        collfreq_ei(ix)
                = coeff * collfreq_elec * density(iion(), ix) / density(ielec(), ix)
                  * (1. + mass_ratio)
                  / std::
                          pow(1. + mass_ratio * temperature(iion(), ix) / temperature(ielec(), ix),
                              1.5);
    });
}

/**
 * Computes the momentum and energy exchange terms between ions and electrons
 */
void compute_momentum_energy_exchange(
        DSpanX momentum_exchange_ei,
        DSpanX momentum_exchange_ie,
        DSpanX energy_exchange_ei,
        DSpanX energy_exchange_ie,
        DViewX collfreq_ei,
        DViewSpX density,
        DViewSpX mean_velocity,
        DViewSpX temperature)
{
    double const mass_ratio(mass(ielec()) / mass(iion()));
    double const me_on_memi(mass(ielec()) / (mass(ielec()) + mass(iion())));
    for_each(policies::parallel_host, collfreq_ei.domain(), [&](IndexX const ix) {
        // momentum exchange terms
        momentum_exchange_ei(ix) = -collfreq_ei(ix) * density(ielec(), ix)
                                   * (mean_velocity(ielec(), ix)
                                      - std::sqrt(mass_ratio) * mean_velocity(iion(), ix));
        momentum_exchange_ie(ix) = -std::sqrt(mass_ratio) * momentum_exchange_ei(ix);

        // energy exchange terms
        energy_exchange_ei(ix) = -3. * collfreq_ei(ix) * me_on_memi * density(ielec(), ix)
                                         * (temperature(ielec(), ix) - temperature(iion(), ix))
                                 - mean_velocity(ielec(), ix) * momentum_exchange_ei(ix);
        energy_exchange_ie(ix)
                = -energy_exchange_ei(ix)
                  - (mean_velocity(ielec(), ix) - std::sqrt(mass_ratio) * mean_velocity(iion(), ix))
                            * momentum_exchange_ei(ix);
    });
}
