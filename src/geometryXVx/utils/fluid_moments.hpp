// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>
#include <quadrature.hpp>

/**
 * Computes fluid moments of the distribution function 
 * Density, mean velocity and temperature.
 */
class FluidMoments
{
private:
    Quadrature<IDimVx> m_integrate_v;

public:
    struct MomentDensity
    {
    };

    struct MomentVelocity
    {
    };

    struct MomentTemperature
    {
    };

    static constexpr MomentDensity s_density = MomentDensity();
    static constexpr MomentVelocity s_velocity = MomentVelocity();
    static constexpr MomentTemperature s_temperature = MomentTemperature();

    FluidMoments(Quadrature<IDimVx> integrate_v);

    ~FluidMoments() = default;

    void operator()(DSpanSpX density, DViewSpXVx allfdistribu, MomentDensity);

    void operator()(
            DSpanSpX mean_velocity,
            DViewSpXVx allfdistribu,
            DViewSpX density,
            MomentVelocity);

    void operator()(
            DSpanSpX temperature,
            DViewSpXVx allfdistribu,
            DViewSpX density,
            DViewSpX mean_velocity,
            MomentTemperature);
};
