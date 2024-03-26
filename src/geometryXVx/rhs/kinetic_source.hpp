#pragma once

#include <cmath>

#include <geometry.hpp>

#include "irighthandside.hpp"

/**
 * @brief A class that describes a source of particles.
 * 
 * The KineticSource class solves the following evolution equation: 
 * df/dt = S_kin
 * Where S_kin = spatial_extent(x) * velocity_shape(v)
 * Since S_kin does not depend on time, we have f(t+dt) = f(t) + S_kin*dt
 * as a solution of this evolution equation. 
 * 
 * spatial_extent defines the location where the source is active.
 * spatial_extent is normalized, so that its integral along the 
 * spatial direction equals one. It has a hyperbolic tangent shape. 
 * It is equal to one in a central zone of the plasma of width defined by the extent parameter.
*  
   velocity_shape defines the velocity profile of the source 
 * in the parallel velocity direction. It is the sum of a source that 
 * injects only density, and a source that injects only energy. If the density 
 * and energy parameters are equal to one (usual case), the resulting velocity_shape 
 * is maxwellian. 
 * 
 * The complete description of the operator can be found in [rhs docs](https://github.com/gyselax/gyselalibxx/blob/main/doc/geometryXVx/kinetic_source.pdf). 
 */
class KineticSource : public IRightHandSide
{
private:
    double m_amplitude;
    double m_density;
    double m_energy;
    double m_temperature;
    host_t<DFieldX> m_spatial_extent;
    host_t<DFieldVx> m_velocity_shape;

public:
    /**
     * @brief Creates an instance of the KineticSource class.
     * @param[in] gridx The mesh in the x direction. 
     * @param[in] gridv The mesh in the vx direction. 
     * @param[in] extent A parameter that sets the spatial extent of the source. 
     * @param[in] stiffness A parameter that sets stiffness of the source extent. 
     * @param[in] amplitude A parameter that sets the amplitude of the source. 
     * @param[in] density A parameter that sets the density of the source. 
     * @param[in] energy A parameter that sets the energy of the source. 
     * @param[in] temperature A parameter that sets the temperature of the source. 
     */
    KineticSource(
            IDomainX const& gridx,
            IDomainVx const& gridv,
            double extent,
            double stiffness,
            double amplitude,
            double density,
            double energy,
            double temperature);

    ~KineticSource() override = default;

    /**
     * @brief Update the distribution function following the KineticSource operator. 
     *
     * Update the distribution function for both electrons and ions to show how
     * it is modified following the effect of the KineticSource operator.
     *
     * @param[inout] allfdistribu The distribution function.
     * @param[in] dt The time step over which the collisions occur.
     *
     * @return A span referencing the distribution function passed as argument.
     */
    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;
};
