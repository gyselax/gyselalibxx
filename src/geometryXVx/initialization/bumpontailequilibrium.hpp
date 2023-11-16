// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>
#include <species_info.hpp>

#include "iequilibrium.hpp"

/// Equilibrium operator as the sum of two Maxwellian. This initializes all species.
class BumpontailEquilibrium : public IEquilibrium
{
    // density of the bump-on-tail part for all kinetic species
    FieldSp<double> m_epsilon_bot;

    // temperature of the bump-on-tail for all kinetic species
    FieldSp<double> m_temperature_bot;

    // mean velocity of the bump-on-tail for all kinetic species
    FieldSp<double> m_mean_velocity_bot;

public:
    /*
      Computation of the Maxwellian fM which is the equilibrium part 
      of the distribution function : 
        fM(v) = f1(v) + f2(v) 
      with 
        f1(v) = (1-epsilon)/(sqrt(2*PI))*exp(-v**2/2)  
        f2(v) = epsilon/sqrt(2*PI*T0)*0.5*[exp(-(v-v0)**2/2*T0)
                                           +exp(-(v+v0)**2/2*T0)]
     */
    void compute_twomaxwellian(
            device_t<DSpanVx> fMaxwellian,
            double epsilon_bot,
            double temperature_bot,
            double mean_velocity_bot) const;

    BumpontailEquilibrium(DViewSp epsilon_bot, DViewSp temperature_bot, DViewSp mean_velocity_bot);

    ~BumpontailEquilibrium() override = default;

    device_t<DSpanSpVx> operator()(device_t<DSpanSpVx> allfequilibrium) const override;

    ViewSp<double> epsilon_bot() const
    {
        return m_epsilon_bot;
    }

    ViewSp<double> temperature_bot() const
    {
        return m_temperature_bot;
    }

    ViewSp<double> mean_velocity_bot() const
    {
        return m_mean_velocity_bot;
    }
};
