// SPDX-License-Identifier: MIT
#pragma once

#include "geometry.hpp"
#include "irighthandside.hpp"

/**
 * @brief A class that describes a source of particles.
 * 
 * The KrookSourceConstant class solves the following evolution equation: 
 * df/dt = -amplitude * mask * (f - ftarget)
 * 
 * mask defines the spatial region where the operator is active.
 * 
 * ftarget is a maxwellian characterized by density and temperature, and a zero fluid velocity.
 *
 * amplitude is a constant for both species. 
 * 
 * The solution of the evolution equation is therefore : 
 * f(t+dt) = ftarget + (f(t)-ftarget)*exp(-amplitude*mask*dt)
 */
class KrookSourceConstant : public IRightHandSide
{
private:
    RhsType m_type;
    double m_extent;
    double m_stiffness;
    double m_amplitude;
    double m_density;
    double m_temperature;
    DFieldMemX m_mask;
    DFieldMemVx m_ftarget;

public:
    /**
     * @brief Creates an instance of the KrookSourceConstant class.
     * @param[in] gridx The mesh in the x direction. 
     * @param[in] gridv The mesh in the vx direction. 
     * @param[in] type A RhsType parameter that defines the region where the operator is active. 
     *                 If type = Source, the mask equals one in the central zone of the plasma of width extent; 
                       If type = Sink, the mask equals zero in the central zone of the plasma of width extent;
     * @param[in] extent A parameter that sets the extent of the source. 
     * @param[in] stiffness A parameter that sets the stiffness of the source extent. 
     * @param[in] amplitude A parameter that sets the the amplitude of the source. 
     * @param[in] density A parameter that sets the density of the Maxwellian ftarget. 
     * @param[in] temperature A parameter that sets the temperature of the Maxwellian ftarget. 
     */
    KrookSourceConstant(
            IdxRangeX const& gridx,
            IdxRangeVx const& gridv,
            RhsType const type,
            double extent,
            double stiffness,
            double amplitude,
            double density,
            double temperature);

    /**
     * @brief Creates an instance of the KrookSourceConstant class.
     */
    KrookSourceConstant(KrookSourceConstant&&) = default;

    ~KrookSourceConstant() override = default;

    /**
     * @brief Update the distribution function following the KrookSourceConstant operator. 
     *
     * Update the distribution function for both electrons and ions to show how
     * it is modified following the effect of the KrookSourceConstant operator.
     *
     * @param[inout] allfdistribu The distribution function.
     * @param[in] dt The time step.
     *
     * @return A field referencing the distribution function passed as argument.
     */
    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, double dt) const override;
};
