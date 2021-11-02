#pragma once

#include "geometry.hpp"

/*
 Computing the non-centered Maxwellian function as
   fM(v) = n/(sqrt(2*PI*T))*exp(-(v-u)**2/(2*T))
  with n the density and T the temperature and
  where u is the mean velocity
*/
void maxwellian_initialization(
        double const density,
        double const temperature,
        double const mean_velocity,
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
            IDomainSpXVx const& domSpXVx)
        : m_charge(std::move(charge))
        , m_mass(std::move(mass))
        , m_density_eq(std::move(n_eq))
        , m_temperature_eq(std::move(T_eq))
        , m_mean_velocity_eq(std::move(u_eq))
        , m_maxw_values(select<IDimSp, IDimVx>(domSpXVx))
    {
        for (IndexSp isp : get_domain<IDimSp>(charge)) {
            // Initialization of the Maxwellian --> fill m_maxw_values
            maxwellian_initialization(
                    m_density_eq(isp),
                    m_temperature_eq(isp),
                    m_mean_velocity_eq(isp),
                    m_maxw_values[isp]);
        }
    }

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
