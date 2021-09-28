#pragma once

#include "geometry.h"

class DistributionFunction
{
public:
    // charge of the particles
    int const charge;

    // mass of the particles
    double const mass;

    // equilibrium density
    double const density_eq;

    // equilibrium temperature
    double const temperature_eq;

    // equilibrium mean velocity
    double const mean_velocity_eq;

    // Maxwellian values
    DBlockVx Maxw_values;

    // values of the function
    DBlockXVx values;

    /*
    /// number of particles
    double nb_particles() const;

    /// Kinetic energy
    double kinetic_energy() const;

    /// Entropy
    double entropy() const;

    /// L1-norm
    double l1_norm() const;

    /// L2-norm
    double l2_norm() const;

    /// density n(x)
    DBlockX density() const;

    /// velocity u(x)
    DBlockX velocity() const;

    /// temperature T(x)
    DBlockX temperature() const;

    /// Stress stress(x)
    DBlockX stress() const;

    //Boundary_conditions
    */

public:
    DistributionFunction(
            int const species_charge,
            double const species_mass,
            double const species_n_eq,
            double const species_T_eq,
            double const species_u_eq,
            MDomainXVx const& domXVx);

    void init();

    MDomainXVx domain() const
    {
        return values.domain();
    }

    DViewVx maxw_values() const
    {
        return Maxw_values;
    }
};
