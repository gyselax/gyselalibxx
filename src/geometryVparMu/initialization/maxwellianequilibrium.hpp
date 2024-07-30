// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>
#include <species_info.hpp>

#include "geometry.hpp"
#include "iequilibrium.hpp"
#include "paraconfpp.hpp"

/// Equilibrium operator as Maxwellian. This initializes all species.
class MaxwellianEquilibrium : public IEquilibrium
{
    // mass of all kinetic species
    host_t<DFieldSp> m_mass;

    // equilibrium density of all kinetic species
    host_t<DFieldSp> m_density_eq;

    // equilibrium temperature of all kinetic species
    host_t<DFieldSp> m_temperature_eq;

    // equilibrium mean velocity of all kinetic species
    host_t<DFieldSp> m_mean_velocity_eq;

private:
    // magnetic field
    double m_magnetic_field;

public:
    /**
     * @brief The constructor for the MaxwellianEquilibrium class.
     * @param[in] mass The mass of the species
     * @param[in] density_eq The density of the Maxwellian
     * @param[in] temperature_eq The temperature of the Maxwellian
     * @param[in] mean_velocity_eq The mean velocity of the Maxwellian
     * @param[in] magnetic_field The magnetic field
     */
    MaxwellianEquilibrium(
            host_t<DFieldSp> mass,
            host_t<DFieldSp> density_eq,
            host_t<DFieldSp> temperature_eq,
            host_t<DFieldSp> mean_velocity_eq,
            double magnetic_field);

    ~MaxwellianEquilibrium() override = default;

    /**
     * @brief Read the density, temperature and mean velocity required to initialize the Maxwellian in a YAML input file.
     * @param[in] dom_kinsp Discrete Domain for the kinetic species
     * @param[in] yaml_input_file YAML input file
     * @return an instance of Maxwellian distribution function.
     */
    static MaxwellianEquilibrium init_from_input(
            IDomainSp dom_kinsp,
            PC_tree_t const& yaml_input_file);

    /**
     * @brief Initializes allfequilibrium as a Maxwellian.
     * @param[out] allfequilibrium A Span containing a Maxwellian distribution function.
     * @return A Span containing a Maxwellian distribution function.
     */
    DSpanSpVparMu operator()(DSpanSpVparMu allfequilibrium) const override;


    /**
     * @brief Compute a Maxwellian distribution function.
     * The Maxwellian distribution function is defined as 
     * Compute $fM(v,mu) = (2*PI*T)**1.5*n*exp(-E)$ with
     *  - $n$ the density, $T$ the temperature and $u$ the mean velocity
     *  - $B$ the magnetic field and
     *  - $E$ the energy defined as $E = (0.5*(v-u)**2+mu*B)/T$.
     * @param[out] fMaxwellian A Maxwellian distribution function. 
     * @param[in] mass Mass of the species.
     * @param[in] density A parameter that represents the density of Maxwellian. 
     * @param[in] temperature A parameter that represents the temperature of Maxwellian. 
     * @param[in] mean_velocity A parameter that represents the mean velocity of Maxwellian. 
     * @param[in] magnetic_field Magnetic field.
     */
    static void compute_maxwellian(
            DSpanVparMu const fMaxwellian,
            double const mass,
            double const density,
            double const temperature,
            double const mean_velocity,
            double const magnetic_field);

    /**
     * @brief A method for accessing the m_mass member variable of the class.
     * @return A view containing the m_mass value. 
     */
    host_t<DViewSp> mass() const
    {
        return m_mass;
    }

    /**
     * @brief A method for accessing the m_density_eq member variable of the class.
     * @return A view containing the m_density_eq value. 
     */
    host_t<DViewSp> density_eq() const
    {
        return m_density_eq;
    }

    /**
     * @brief A method for accessing the m_temperature_eq member variable of the class.
     * @return A view containing the m_temperature_eq value. 
     */
    host_t<DViewSp> temperature_eq() const
    {
        return m_temperature_eq;
    }

    /**
     * @brief A method for accessing the m_mean_velocity_eq member variable of the class.
     * @return A view containing the m_velocity_eq value. 
     */
    host_t<DViewSp> mean_velocity_eq() const
    {
        return m_mean_velocity_eq;
    }
};
