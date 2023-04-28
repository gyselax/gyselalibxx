
#include <cassert>
#include <cmath>
#include <optional>

#include <ddc_helper.hpp>
#include <quadrature.hpp>
#include <species_info.hpp>
#include <trapezoid_quadrature.hpp>

#include "collisions_utils.hpp"

namespace {

/**
 * Useful function for computing collision related quantities:
 *  - kernel maxwellian moments
 * Warning: only meaningful for the collision operator!
 */
IndexSp find_ion(IDomainSp const dom_sp)
{
    assert(dom_sp.size() == 2);
    std::optional<IndexSp> iion_opt;
    for (IndexSp const isp : dom_sp) {
        if (charge(isp) > 0) {
            iion_opt = isp;
        }
    }
    return iion_opt.value();
}

} // namespace

/**
 * Computes the spatial profile of nustar, which is constant here,
 * but could be space dependent (for instance to have no collisions in some specific
 * parts of the simulation box).
 */
void compute_nustar_profile(DSpanSpX nustar_profile, double nustar0)
{
    double const Lx = ddcHelper::total_interval_length(ddc::get_domain<IDimX>(nustar_profile));
    ddc::for_each(ddc::policies::parallel_host, nustar_profile.domain(), [&](IndexSpX const ispx) {
        double const coeff = std::sqrt(mass(ielec()) / mass(ddc::select<IDimSp>(ispx)))
                             * std::pow(charge(ddc::select<IDimSp>(ispx)), 4) / Lx;
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
    ddc::for_each(ddc::policies::parallel_host, collfreq.domain(), [&](IndexSpX const ispx) {
        collfreq(ispx) = nustar_profile(ispx) * density(ispx) / std::pow(temperature(ispx), 1.5);
    });
}

/**
 * Computes the two species collision frequency collfreq_ei
 */
void compute_collfreq_ab(
        DSpanSp collfreq_ab,
        DViewSp nustar_profile,
        DViewSp density,
        DViewSp temperature)
{
    IndexSp const iion = find_ion(density.domain<IDimSp>());
    double const charge_ratio(charge(iion) / charge(ielec()));
    double const me_on_mi(mass(ielec()) / mass(iion));
    double const collfreq_elec(
            nustar_profile(ielec()) * density(ielec()) / std::pow(temperature(ielec()), 1.5));

    collfreq_ab(ielec())
            = std::sqrt(2.) * charge_ratio * charge_ratio * collfreq_elec * density(iion)
              / density(ielec()) * (1. + me_on_mi)
              / std::pow(1. + me_on_mi * temperature(iion) / temperature(ielec()), 1.5);

    double const collfreq_ion(
            nustar_profile(iion) * density(iion) / std::pow(temperature(iion), 1.5));
    collfreq_ab(iion)
            = std::sqrt(2.) / (charge_ratio * charge_ratio) * collfreq_ion * density(ielec())
              / density(iion) * (1. + 1. / me_on_mi)
              / std::pow(1. + 1. / me_on_mi * temperature(ielec()) / temperature(iion), 1.5);
}

/**
 * Computes the momentum and energy exchange terms between ions and electrons
 */
void compute_momentum_energy_exchange(
        DSpanSp momentum_exchange_ab,
        DSpanSp energy_exchange_ab,
        DViewSp collfreq_ab,
        DViewSp density,
        DViewSp mean_velocity,
        DViewSp temperature)
{
    IndexSp const iion = find_ion(density.domain<IDimSp>());
    double const mass_ratio(mass(ielec()) / mass(iion));
    double const me_on_memi(mass(ielec()) / (mass(ielec()) + mass(iion)));
    // momentum exchange terms
    momentum_exchange_ab(ielec())
            = -collfreq_ab(ielec()) * density(ielec())
              * (mean_velocity(ielec()) - std::sqrt(mass_ratio) * mean_velocity(iion));
    momentum_exchange_ab(iion) = -std::sqrt(mass_ratio) * momentum_exchange_ab(ielec());

    // energy exchange terms
    energy_exchange_ab(ielec()) = -3. * collfreq_ab(ielec()) * me_on_memi * density(ielec())
                                          * (temperature(ielec()) - temperature(iion))
                                  - mean_velocity(ielec()) * momentum_exchange_ab(ielec());
    energy_exchange_ab(iion)
            = -energy_exchange_ab(ielec())
              - (mean_velocity(ielec()) - std::sqrt(mass_ratio) * mean_velocity(iion))
                        * momentum_exchange_ab(ielec());
}
