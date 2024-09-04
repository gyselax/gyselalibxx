// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "iequilibrium.hpp"
#include "species_info.hpp"

/// Equilibrium operator as normalized Maxwellian. This initializes all species.
class MaxwellianEquilibrium : public IEquilibrium
{
    // equilibrium density of all kinetic species
    DFieldMemSpTor2D m_density_eq;

    // equilibrium temperature of all kinetic species
    DFieldMemSpTor2D m_temperature_eq;

    // equilibrium mean velocity of all kinetic species
    DFieldMemSpTor2D m_mean_velocity_eq;

private:
    // magnetic field
    DFieldMemTor2D m_magnetic_field;

public:
    /**
     * @brief The constructor for the MaxwellianEquilibrium class.
     * @param[in] density_eq The density of the Maxwellian
     * @param[in] temperature_eq The temperature of the Maxwellian
     * @param[in] mean_velocity_eq The mean velocity of the Maxwellian
     * @param[in] magnetic_field The magnetic field
     */
    MaxwellianEquilibrium(
            host_t<DFieldSpTor2D> density_eq,
            host_t<DFieldSpTor2D> temperature_eq,
            host_t<DFieldSpTor2D> mean_velocity_eq,
            host_t<DFieldTor2D> magnetic_field);

    ~MaxwellianEquilibrium() override = default;

    /**
     * @brief Initializes allfequilibrium as a Maxwellian.
     * @param[out] allfequilibrium A Field containing a Maxwellian distribution function.
     * @return A Field containing a Maxwellian distribution function.
     */
    DFieldSpV2DTor2D operator()(DFieldSpV2DTor2D allfequilibrium) const override;


    /**
     * @brief Compute a Maxwellian distribution function.
     * The normalized Maxwellian distribution function is defined as 
     *  $fM(r,theta,vpar,mu) = (2*PI*T(r,theta))**(-1.5)*n(r,theta)*exp(-E)$ with
     *  - $n$ the density, $T$ the temperature and $u$ the mean velocity
     *  - $B$ the magnetic field and
     *  - $E$ the energy defined as $E = (0.5*(v-u)**2+mu*B)/T$.
     * @param[out] fMaxwellian A Maxwellian distribution function. 
     * @param[in] density A parameter that represents the density of Maxwellian. 
     * @param[in] temperature A parameter that represents the temperature of Maxwellian. 
     * @param[in] mean_velocity A parameter that represents the mean velocity of Maxwellian. 
     * @param[in] magnetic_field Magnetic field.
     */
    static void compute_maxwellian(
            DFieldV2DTor2D const fMaxwellian,
            DConstFieldTor2D density,
            DConstFieldTor2D temperature,
            DConstFieldTor2D mean_velocity,
            DConstFieldTor2D magnetic_field);

    /**
     * @brief A method for accessing the m_density_eq member variable of the class.
     * @return A view containing the m_density_eq value. 
     */
    DConstFieldSpTor2D density_eq() const
    {
        return m_density_eq;
    }

    /**
     * @brief A method for accessing the m_temperature_eq member variable of the class.
     * @return A view containing the m_temperature_eq value. 
     */
    DConstFieldSpTor2D temperature_eq() const
    {
        return m_temperature_eq;
    }

    /**
     * @brief A method for accessing the m_mean_velocity_eq member variable of the class.
     * @return A view containing the m_velocity_eq value. 
     */
    DConstFieldSpTor2D mean_velocity_eq() const
    {
        return m_mean_velocity_eq;
    }
};
