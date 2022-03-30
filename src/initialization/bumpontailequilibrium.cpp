// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "bumpontailequilibrium.hpp"

using std::sqrt, std::exp;

BumpontailEquilibrium::BumpontailEquilibrium(
        DViewSp const epsilon_bot,
        DViewSp const temperature_bot,
        DViewSp const mean_velocity_bot)
    : m_epsilon_bot(epsilon_bot)
    , m_temperature_bot(temperature_bot)
    , m_mean_velocity_bot(mean_velocity_bot)
{
}

DSpanSpVx BumpontailEquilibrium::operator()(DSpanSpVx const allfequilibrium) const
{
    IDomainVx const gridvx = allfequilibrium.domain<IDimVx>();
    IDomainSp const gridsp = allfequilibrium.domain<IDimSp>();

    // Initialization of the maxwellian
    DFieldVx maxwellian(gridvx);
    for_each(gridsp, [&](IndexSp const isp) {
        compute_twomaxwellian(
                maxwellian,
                m_epsilon_bot(isp),
                m_temperature_bot(isp),
                m_mean_velocity_bot(isp));

        for_each(gridvx, [&](IndexVx const iv) { allfequilibrium(isp, iv) = maxwellian(iv); });
    });
    return allfequilibrium;
}

/*
Computation of the Maxwellian fM which is the equilibrium part 
of the distribution function : 
  fM(v) = f1(v) + f2(v) 
with 
  f1(v) = (1-epsilon)/(sqrt(2*PI))*exp(-v**2/2)  
  f2(v) = epsilon/sqrt(2*PI*T0)[exp(-(v-v0)**2/2*T0)]
*/
void BumpontailEquilibrium::compute_twomaxwellian(
        DSpanVx const fMaxwellian,
        double const epsilon_bot,
        double const temperature_bot,
        double const mean_velocity_bot) const
{
    double const inv_sqrt_2pi = 1. / sqrt(2. * M_PI);
    double const norm_f2 = inv_sqrt_2pi / sqrt(temperature_bot);
    IDomainVx const gridvx = fMaxwellian.domain();
    for (IndexVx const iv : gridvx) {
        CoordVx const v = to_real(iv);
        // bulk plasma particles
        double const f1_v = (1. - epsilon_bot) * inv_sqrt_2pi * exp(-0.5 * v * v);
        // beam
        double const f2_v = epsilon_bot * norm_f2
                            * (exp(-(v - mean_velocity_bot) * (v - mean_velocity_bot)
                                   / (2. * temperature_bot)));
        // fM(v) = f1(v) + f2(v)
        fMaxwellian(iv) = f1_v + f2_v;
    }
}
