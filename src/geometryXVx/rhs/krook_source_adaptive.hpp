#pragma once

#include <geometry.hpp>

#include "irighthandside.hpp"

/**
 * @brief A class that describes a source of particles.
 * 
 * The KrookSourceAdaptive class solves the following evolution equation: 
 * df/dt = -amplitude * mask * (f - ftarget)
 * 
 * mask defines the spatial region where the operator is active.
 * 
 * ftarget is a maxwellian characterized by density and temperature, and a zero fluid velocity.
 *
 * amplitude depends on space, time and the considered species so that: 
 * - amplitude(ions) = m_amplitude = constant  
 * - amplitude(electrons, x, t) = m_amplitude (density_ions(x,t) - m_density) / (density_electrons(x,t) - m_density)  
 * so that the operator conserves locally the charge. 
 * 
 * The complete description of the operator can be found in [rhs docs](https://github.com/gyselax/gyselalibxx/blob/main/src/geometryXVx/rhs/doc/krook_source.pdf). 
 */
class KrookSourceAdaptive : public IRightHandSide
{
private:
    RhsType m_type;
    double m_extent;
    double m_stiffness;
    double m_amplitude;
    double m_density;
    double m_temperature;
    DFieldX m_mask;
    DFieldVx m_ftarget;

public:
    /**
     * @brief Creates an instance of the KrookSourceAdaptive class.
     * @param[in] gridx The mesh in the x direction. 
     * @param[in] gridvx The mesh in the vx direction. 
     * @param[in] type A RhsType parameter that defines the region where the operator is active. 
     *                 If type = Source, the mask equals one in the central zone of the plasma of width extent; 
                       If type = Sink, the mask equals zero in the central zone of the plasma of width extent;
     * @param[in] extent A parameter that sets the extent of the source. 
     * @param[in] stiffness A parameter that sets the stiffness of the source extent. 
     * @param[in] amplitude A parameter that sets the amplitude of the source. 
     * @param[in] density A parameter that sets the density of the Maxwellian ftarget. 
     * @param[in] temperature A parameter that sets the temperature of the Maxwellian ftarget. 
     */
    KrookSourceAdaptive(
            IDomainX const& gridx,
            IDomainVx const& gridvx,
            RhsType const type,
            double extent,
            double stiffness,
            double amplitude,
            double density,
            double temperature);

    /**
     * @brief Creates an instance of the KrookSourceAdaptive class.
     */
    KrookSourceAdaptive(KrookSourceAdaptive&&) = default;

    ~KrookSourceAdaptive() override = default;

    /**
     * @brief Update the distribution function following the KrookSourceAdaptive operator. 
     *
     * Update the distribution function for both electrons and ions to show how
     * it is modified following the effect of the KrookSourceAdaptive operator.
     *
     * @param[inout] allfdistribu The distribution function.
     * @param[in] dt The time step.
     *
     * @return A span referencing the distribution function passed as argument.
     */
    device_t<DSpanSpXVx> operator()(device_t<DSpanSpXVx> allfdistribu, double dt) const override;

public:
    /**
     * @brief Computes the amplitude coefficient of the KrookSourceAdaptive operator. 
     *
     * This coefficient depends on the considered species and ensures that 
     * the operator conserves the charge locally. 
     *
     * @param[inout] amplitudes A Span that contains on output the amplitude 
     *                         coefficients for each species. 
     * @param[in] allfdistribu The distribution function.
     */
    void get_amplitudes(device_t<DSpanSpX> amplitudes, device_t<DViewSpXVx> allfdistribu) const;

    /**
     * @brief Computes the expression of the time derivative of the distribution function. 
     * 
     * The expression is df = -amplitude * mask * (f - ftarget). This function
     * is used for the time integrator (RK2 for instance).
     *
     * @param[inout] df The time derivative.
     * @param[in] f The distribution function.
     * @param[in] f0 An optional parameter.
     */
    void get_derivative(DSpanSpXVx df, DViewSpXVx f, DViewSpXVx f0) const;
};
