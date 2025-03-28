

# File kinetic\_source.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**kinetic\_source.hpp**](kinetic__source_8hpp.md)

[Go to the documentation of this file](kinetic__source_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>

#include "geometry.hpp"
#include "irighthandside.hpp"

class KineticSource : public IRightHandSide
{
private:
    double m_amplitude;
    double m_density;
    double m_energy;
    double m_temperature;
    DFieldMemX m_spatial_extent;
    DFieldMemVx m_velocity_shape;

public:
    KineticSource(
            IdxRangeX const& gridx,
            IdxRangeVx const& gridv,
            double extent,
            double stiffness,
            double amplitude,
            double density,
            double energy,
            double temperature);

    ~KineticSource() override = default;

    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, double dt) const override;
};
```


