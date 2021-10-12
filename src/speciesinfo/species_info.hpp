#pragma once

#include "geometry.h"

class SpeciesInformation
{
private:
    // charge of the particles
    BlockSp<int> const m_charge;

    // mass of the particles
    BlockSp<double> const m_mass;

    // equilibrium density
    BlockSp<double> const m_density_eq;

    // equilibrium temperature
    BlockSp<double> const m_temperature_eq;

    // equilibrium mean velocity
    BlockSp<double> const m_mean_velocity_eq;

    // Maxwellian values
    BlockSpVx<double> const m_maxw_values;

    // initial perturbed mode
    BlockSp<int> const m_init_perturb_mode;

    // initial perturbation amplitude
    BlockSp<double> const m_init_perturb_amplitude;

public:
    SpeciesInformation(
            BlockSp<int> charge,
            BlockSp<double> mass,
            BlockSp<double> n_eq,
            BlockSp<double> T_eq,
            BlockSp<double> u_eq,
            BlockSp<int> init_perturb_mode,
            BlockSp<double> init_perturb_amplitude,
            MDomainSpXVx const& domSpXVx)
        : m_charge(std::move(charge))
        , m_mass(std::move(mass))
        , m_density_eq(std::move(n_eq))
        , m_temperature_eq(std::move(T_eq))
        , m_mean_velocity_eq(std::move(u_eq))
        , m_maxw_values(select<MeshSp, MeshVx>(domSpXVx))
        , m_init_perturb_mode(std::move(init_perturb_mode))
        , m_init_perturb_amplitude(std::move(init_perturb_amplitude))
    {
    }

    MCoordSp ielec() const
    {
        return MCoordSp(0);
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

    void init(DSpanSpXVx fdistribu);
};
