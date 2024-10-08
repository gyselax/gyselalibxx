
#include <cassert>
#include <cmath>
#include <optional>

#include "collisions_utils.hpp"
#include "ddc_helper.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"
#include "trapezoid_quadrature.hpp"

namespace {

/**
 * Useful function for computing collision related quantities:
 *  - kernel maxwellian moments
 * Warning: only meaningful for the collision operator!
 */
IdxSp find_ion(IdxRangeSp const idx_range_sp)
{
    assert(idx_range_sp.size() == 2);
    std::optional<IdxSp> iion_opt;
    for (IdxSp const isp : idx_range_sp) {
        if (charge(isp) > 0.) {
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
void compute_nustar_profile(DFieldSpX nustar_profile, double nustar0)
{
    double const Lx = ddcHelper::total_interval_length(get_idx_range<GridX>(nustar_profile));

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(nustar_profile),
            KOKKOS_LAMBDA(IdxSpX const ispx) {
                double const coeff = Kokkos::sqrt(mass(ielec()) / mass(ddc::select<Species>(ispx)))
                                     * Kokkos::pow(charge(ddc::select<Species>(ispx)), 4) / Lx;
                nustar_profile(ispx) = coeff * nustar0;
            });
}

/**
 * Computes the space and species dependent collision frequency collfreq.
 */
void compute_collfreq(
        DFieldSpX collfreq,
        DConstFieldSpX nustar_profile,
        DConstFieldSpX density,
        DConstFieldSpX temperature)
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(collfreq),
            KOKKOS_LAMBDA(IdxSpX const ispx) {
                collfreq(ispx) = nustar_profile(ispx) * density(ispx)
                                 / Kokkos::pow(temperature(ispx), 1.5);
            });
}

/**
 * Computes the two species collision frequency collfreq_ei
 */
void compute_collfreq_ab(
        DFieldSpX collfreq_ab,
        DConstFieldSpX nustar_profile,
        DConstFieldSpX density,
        DConstFieldSpX temperature)
{
    IdxSp const iion = find_ion(get_idx_range<Species>(density));
    double const charge_ratio(charge(iion) / charge(ielec()));
    double const me_on_mi(mass(ielec()) / mass(iion));

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range<GridX>(collfreq_ab),
            KOKKOS_LAMBDA(IdxX const ix) {
                double const collfreq_elec(
                        nustar_profile(ielec(), ix) * density(ielec(), ix)
                        / Kokkos::pow(temperature(ielec(), ix), 1.5));

                collfreq_ab(ielec(), ix) = Kokkos::sqrt(2.) * charge_ratio * charge_ratio
                                           * collfreq_elec * density(iion, ix)
                                           / density(ielec(), ix) * (1. + me_on_mi)
                                           / Kokkos::
                                                   pow(1.
                                                               + me_on_mi * temperature(iion, ix)
                                                                         / temperature(ielec(), ix),
                                                       1.5);

                double const collfreq_ion(
                        nustar_profile(iion, ix) * density(iion, ix)
                        / Kokkos::pow(temperature(iion, ix), 1.5));
                collfreq_ab(iion, ix)
                        = Kokkos::sqrt(2.) / (charge_ratio * charge_ratio) * collfreq_ion
                          * density(ielec(), ix) / density(iion, ix) * (1. + 1. / me_on_mi)
                          / Kokkos::
                                  pow(1.
                                              + 1. / me_on_mi * temperature(ielec(), ix)
                                                        / temperature(iion, ix),
                                      1.5);
            });
}

/**
 * Computes the momentum and energy exchange terms between ions and electrons
 */
void compute_momentum_energy_exchange(
        DFieldSpX momentum_exchange_ab,
        DFieldSpX energy_exchange_ab,
        DConstFieldSpX collfreq_ab,
        DConstFieldSpX density,
        DConstFieldSpX mean_velocity,
        DConstFieldSpX temperature)
{
    IdxSp const iion = find_ion(get_idx_range<Species>(density));
    double const mass_ratio(mass(ielec()) / mass(iion));
    double const me_on_memi(mass(ielec()) / (mass(ielec()) + mass(iion)));
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range<GridX>(collfreq_ab),
            KOKKOS_LAMBDA(IdxX const ix) {
                // momentum exchange terms
                momentum_exchange_ab(ielec(), ix)
                        = -collfreq_ab(ielec(), ix) * density(ielec(), ix)
                          * (mean_velocity(ielec(), ix)
                             - Kokkos::sqrt(mass_ratio) * mean_velocity(iion, ix));
                momentum_exchange_ab(iion, ix)
                        = -Kokkos::sqrt(mass_ratio) * momentum_exchange_ab(ielec(), ix);

                // energy exchange terms
                energy_exchange_ab(ielec(), ix)
                        = -3. * collfreq_ab(ielec(), ix) * me_on_memi * density(ielec(), ix)
                                  * (temperature(ielec(), ix) - temperature(iion, ix))
                          - mean_velocity(ielec(), ix) * momentum_exchange_ab(ielec(), ix);
                energy_exchange_ab(iion, ix)
                        = -energy_exchange_ab(ielec(), ix)
                          - (mean_velocity(ielec(), ix)
                             - Kokkos::sqrt(mass_ratio) * mean_velocity(iion, ix))
                                    * momentum_exchange_ab(ielec(), ix);
            });
}
