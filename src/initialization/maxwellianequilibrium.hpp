// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>
#include <species_info.hpp>

#include "iequilibrium.hpp"

/// Equilibrium operator as Maxwellian. This initializes all species.
class MaxwellianEquilibrium : public IEquilibrium
{
    // equilibrium density of all kinetic species
    FieldSp<double> m_density_eq;

    // equilibrium temperature of all kinetic species
    FieldSp<double> m_temperature_eq;

    // equilibrium mean velocity of all kinetic species
    FieldSp<double> m_mean_velocity_eq;

private:
    /*
 Computing the non-centered Maxwellian function as
   fM(v) = n/(sqrt(2*PI*T))*exp(-(v-u)**2/(2*T))
  with n the density and T the temperature and
  where u is the mean velocity
*/
    void compute_maxwellian(
            DSpanVx fMaxwellian,
            double density,
            double temperature,
            double mean_velocity) const;

public:
    MaxwellianEquilibrium(DViewSp density_eq, DViewSp temperature_eq, DViewSp mean_velocity_eq);

    ~MaxwellianEquilibrium() override = default;

    DSpanSpVx operator()(DSpanSpVx allfequilibrium) const override;

    ViewSp<double> density_eq() const
    {
        return m_density_eq;
    }

    ViewSp<double> temperature_eq() const
    {
        return m_temperature_eq;
    }

    ViewSp<double> mean_velocity_eq() const
    {
        return m_mean_velocity_eq;
    }
};
