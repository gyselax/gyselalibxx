// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>
#include <species_info.hpp>

#include "iequilibrium.hpp"

/**
 * @brief A class that initializes the distribution function as a sum of two Maxwellian functions.
 *
 * This class initializes the distribution function as a sum of two Maxwellian, 
 * enabling the study of the so-called bump-on-tail instability. One of the 
 * Maxwellians represents the bulk of the distribution function that has no mean 
 * velocity, and the other Maxwellian corresponds to high velocity particles. 
 * The second Maxwellian is referred to as the "bump-on-tail" Maxwellian.
 */
class BumpontailEquilibrium : public IEquilibrium
{
    /**Density of the bump-on-tail part for all kinetic species*/
    host_t<DFieldSp> m_epsilon_bot;

    /**Temperature of the bump-on-tail for all kinetic species*/
    host_t<DFieldSp> m_temperature_bot;

    /**Mean velocity of the bump-on-tail for all kinetic species*/
    host_t<DFieldSp> m_mean_velocity_bot;

public:
    /**
     * @brief Compute a distribution function defined as a sum of two Maxwellians.
     * This distribution function can be written as 
     * $f(x,v) = f1(v) + f2(v) $
     * with 
     * $f1(v) = (1-epsilon)/(sqrt(2*PI))*exp(-v**2/2)$
     * $f2(v) = epsilon/sqrt(2*PI*T0)[exp(-(v-v0)**2/2*T0)$
     * @param[out] fMaxwellian The initial distribution function. 
     * @param[in] epsilon_bot A parameter that represents the density of the bump-on-tail Maxwellian. 
     * @param[in] temperature_bot A parameter that represents the temperature of the bump-on-tail Maxwellian. 
     * @param[in] mean_velocity_bot A parameter that represents the mean velocity of the bump-on-tail Maxwellian. 
     */
    void compute_twomaxwellian(
            DSpanVx fMaxwellian,
            double epsilon_bot,
            double temperature_bot,
            double mean_velocity_bot) const;
    /**
     * @brief Creates an instance of the BumpontailEquilibrium class.
     * @param[in] epsilon_bot A parameter that represents the density of the bump-on-tail Maxwellian for each species. 
     * @param[in] temperature_bot A parameter that represents the temperature of the bump-on-tail Maxwellian for each species. 
     * @param[in] mean_velocity_bot A parameter that represents the mean velocity of the bump-on-tail Maxwellian for each species. 
     */
    BumpontailEquilibrium(
            host_t<DFieldSp> epsilon_bot,
            host_t<DFieldSp> temperature_bot,
            host_t<DFieldSp> mean_velocity_bot);

    ~BumpontailEquilibrium() override = default;

    /**
     * @brief Initializes the distribution function as the sum of a bulk and a bump-on-tail Maxwellians. 
     * @param[out] allfequilibrium The initialized distribution function.
     * @return The initialized distribution function.
     */
    DSpanSpVx operator()(DSpanSpVx allfequilibrium) const override;

    /**
     * @brief A method for accessing the m_epsilon_bot member variable of the class.
     * @return a View containing the m_epsilon_bot variable.
     */
    host_t<DViewSp> epsilon_bot() const
    {
        return m_epsilon_bot;
    }

    /**
     * @brief A method for accessing the m_temperature_bot member variable of the class.
     * @return a View containing the m_temperature_bot variable.
     */
    host_t<DViewSp> temperature_bot() const
    {
        return m_temperature_bot;
    }

    /**
     * @brief A method for accessing the m_mean_velocity_bot member variable of the class.
     * @return a View containing the m_velocity_bot variable.
     */
    host_t<DViewSp> mean_velocity_bot() const
    {
        return m_mean_velocity_bot;
    }
};
