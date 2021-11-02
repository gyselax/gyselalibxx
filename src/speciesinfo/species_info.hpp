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
    // charge of the particles
    FieldSp<int> const m_charge;

    // mass of the particles
    FieldSp<double> const m_mass;

    // equilibrium density
    FieldSp<double> const m_density_eq;

    // equilibrium temperature
    FieldSp<double> const m_temperature_eq;

    // equilibrium mean velocity
    FieldSp<double> const m_mean_velocity_eq;

    // Maxwellian values
    FieldSpVx<double> m_maxw_values;

public:
    SpeciesInformation(
            FieldSp<int> charge,
            FieldSp<double> mass,
            FieldSp<double> n_eq,
            FieldSp<double> T_eq,
            FieldSp<double> u_eq,
            IDomainSpXVx const& domSpXVx);

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

    ViewSpVx<double> maxw_values() const
    {
        return m_maxw_values;
    }
};
