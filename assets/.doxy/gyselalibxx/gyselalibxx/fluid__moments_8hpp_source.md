

# File fluid\_moments.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**utils**](dir_8b9ab5da15e50812e4d198d35fde42ae.md) **>** [**fluid\_moments.hpp**](fluid__moments_8hpp.md)

[Go to the documentation of this file](fluid__moments_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "quadrature.hpp"

class FluidMoments
{
private:
    Quadrature<IdxRangeVx, IdxRangeSpXVx> m_integrate_v;

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

    FluidMoments(Quadrature<IdxRangeVx, IdxRangeSpXVx> integrate_v);

    ~FluidMoments() = default;

    void operator()(double& density, DConstFieldVx fdistribu, MomentDensity moment_density);

    void operator()(DFieldSpX density, DConstFieldSpXVx allfdistribu, MomentDensity moment_density);

    void operator()(
            double& mean_velocity,
            DConstFieldVx fdistribu,
            double density,
            MomentVelocity moment_velocity);

    void operator()(
            DFieldSpX mean_velocity,
            DConstFieldSpXVx allfdistribu,
            DConstFieldSpX density,
            MomentVelocity moment_velocity);

    void operator()(
            double& temperature,
            DConstFieldVx fdistribu,
            double density,
            double mean_velocity,
            MomentTemperature moment_temperature);

    void operator()(
            DFieldSpX temperature,
            DConstFieldSpXVx allfdistribu,
            DConstFieldSpX density,
            DConstFieldSpX mean_velocity,
            MomentTemperature moment_temperature);
};
```


