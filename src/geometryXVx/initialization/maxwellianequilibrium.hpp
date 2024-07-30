// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>
#include <species_info.hpp>

#include "geometry.hpp"
#include "iequilibrium.hpp"
#include "paraconfpp.hpp"

/**
 * @brief A class that initializes the distribution function as a Maxwellian.
 */
class MaxwellianEquilibrium : public IEquilibrium
{
    /**Equilibrium density of all kinetic species*/
    host_t<DFieldMemSp> m_density_eq;

    /**Equilibrium temperature of all kinetic species*/
    host_t<DFieldMemSp> m_temperature_eq;

    /**Equilibrium mean velocity of all kinetic species*/
    host_t<DFieldMemSp> m_mean_velocity_eq;

public:
    /**
     * @brief The constructor for the MaxwellianEquilibrium class.
     * @param[in] density_eq The density of the Maxwellian
     * @param[in] temperature_eq The temperature of the Maxwellian
     * @param[in] mean_velocity_eq The mean velocity of the Maxwellian
     */
    MaxwellianEquilibrium(
            host_t<DFieldMemSp> density_eq,
            host_t<DFieldMemSp> temperature_eq,
            host_t<DFieldMemSp> mean_velocity_eq);

    ~MaxwellianEquilibrium() override = default;

    /**
     * @brief Read the density, temperature and mean velocity required to initialize the Maxwellian in a YAML input file.
     * @param[in] dom_kinsp Discrete Domain for the kinetic species
     * @param[in] yaml_input_file YAML input file
     * @return an instance of Maxwellian distribution function.
     */
    static MaxwellianEquilibrium init_from_input(
            IdxRangeSp dom_kinsp,
            PC_tree_t const& yaml_input_file);

    /**
     * @brief Initializes allfequilibrium as a Maxwellian.
     * @param[out] allfequilibrium A Span containing a Maxwellian distribution function.
     * @return A Span containing a Maxwellian distribution function.
     */
    DSpanSpVx operator()(DSpanSpVx allfequilibrium) const override;

    /**
     * @brief Compute a Maxwellian distribution function.
     * The Maxwellian distribution function is defined as 
     * $f_M(v) = n/(sqrt(2*PI*T))*exp(-(v-u)**2/(2*T))$
     * with $n$ the density, $T$ the temperature and
     * $u$ is the mean velocity.
     * @param[out] fMaxwellian A Maxwellian distribution function. 
     * @param[in] density A parameter that represents the density of Maxwellian. 
     * @param[in] temperature A parameter that represents the temperature of Maxwellian. 
     * @param[in] mean_velocity A parameter that represents the mean velocity of Maxwellian. 
     */
    static void compute_maxwellian(
            DSpanVx const fMaxwellian,
            double const density,
            double const temperature,
            double const mean_velocity);

    /**
     * @brief A method for accessing the m_density_eq member variable of the class.
     * @return A view containing the m_density_eq value. 
     */
    host_t<DConstFieldSp> density_eq() const
    {
        return m_density_eq;
    }

    /**
     * @brief A method for accessing the m_temperature_eq member variable of the class.
     * @return A view containing the m_temperature_eq value. 
     */
    host_t<DConstFieldSp> temperature_eq() const
    {
        return m_temperature_eq;
    }

    /**
     * @brief A method for accessing the m_mean_velocity_eq member variable of the class.
     * @return A view containing the m_velocity_eq value. 
     */
    host_t<DConstFieldSp> mean_velocity_eq() const
    {
        return m_mean_velocity_eq;
    }
};
