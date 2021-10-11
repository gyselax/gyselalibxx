#pragma once

#include "geometry.h"

class DistributionFunction
{
public:
    // charge of the particles
    int const m_charge;

    // mass of the particles
    double const m_mass;

    // equilibrium density
    double const m_density_eq;

    // equilibrium temperature
    double const m_temperature_eq;

    // equilibrium mean velocity
    double const m_mean_velocity_eq;

    // Maxwellian values
    DBlockVx m_Maxw_values;

    // initial perturbed mode
    int const m_init_perturb_mode;

    // initial perturbation amplitude
    double const m_init_perturb_amplitude;

    // values of the function
    DBlockXVx m_values;

public:
    DistributionFunction(
            int const charge,
            double const mass,
            double const n_eq,
            double const T_eq,
            double const u_eq,
            int const init_perturb_mode,
            double const init_perturb_amplitude,
            MDomainXVx const& domXVx)
        : m_charge(charge)
        , m_mass(mass)
        , m_density_eq(n_eq)
        , m_temperature_eq(T_eq)
        , m_mean_velocity_eq(u_eq)
        , m_Maxw_values(select<MeshVx>(domXVx))
        , m_init_perturb_mode(init_perturb_mode)
        , m_init_perturb_amplitude(init_perturb_amplitude)
        , m_values(domXVx)
    {
    }

    MDomainX domainX() const
    {
        return m_values.domain<MeshX>();
    }

    MDomainVx domainVx() const
    {
        return m_values.domain<MeshVx>();
    }

    MDomainXVx domain() const
    {
        return m_values.domain();
    }

    DViewVx maxw_values() const
    {
        return m_Maxw_values;
    }

    void init();
};
