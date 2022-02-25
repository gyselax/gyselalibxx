// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/*
 Computing the non-centered Maxwellian function as
   fM(v) = n/(sqrt(2*PI*T))*exp(-(v-u)**2/(2*T))
  with n the density and T the temperature and
  where u is the mean velocity
*/
void maxwellian_initialization(
        double density,
        double temperature,
        double mean_velocity,
        DSpanVx fMaxwellian);

class SpeciesInformation
{
private:
    // charge of the particles (kinetic + adiabatic)
    FieldSp<int> const m_charge;

    // mass of the particles of all kinetic species
    FieldSp<double> const m_mass;

    // equilibrium density of all kinetic species
    FieldSp<double> const m_density_eq;

    // equilibrium temperature of all kinetic species
    FieldSp<double> const m_temperature_eq;

    // equilibrium mean velocity of all kinetic species
    FieldSp<double> const m_mean_velocity_eq;

    // Initial perturbation amplitude of all kinetic species
    FieldSp<double> const m_perturb_amplitude;

    // Initial perturbation mode of all kinetic species
    FieldSp<int> const m_perturb_mode;

    // Maxwellian values of all kinetic species
    FieldSpVx<double> m_maxw_values;

public:
    SpeciesInformation(
            FieldSp<int> charge,
            FieldSp<double> mass,
            FieldSp<double> n_eq,
            FieldSp<double> T_eq,
            FieldSp<double> u_eq,
            FieldSp<double> perturb_amplitude,
            FieldSp<int> perturb_mode,
            IDomainSpXVx const& domSpXVx);

    // Consider that the electron species is always at the 0 position
    IndexSp ielec() const
    {
        return IndexSp(0);
    }

    ViewSp<int> charge() const
    {
        return m_charge;
    }

    ViewSp<double> mass() const
    {
        return m_mass;
    }

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

    ViewSp<double> perturb_amplitude() const
    {
        return m_perturb_amplitude;
    }

    ViewSp<int> perturb_mode() const
    {
        return m_perturb_mode;
    }

    ViewSpVx<double> maxw_values() const
    {
        return m_maxw_values;
    }
};
