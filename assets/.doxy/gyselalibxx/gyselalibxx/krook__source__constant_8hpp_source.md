

# File krook\_source\_constant.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**krook\_source\_constant.hpp**](krook__source__constant_8hpp.md)

[Go to the documentation of this file](krook__source__constant_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "geometry.hpp"
#include "irighthandside.hpp"

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
    KrookSourceConstant(
            IdxRangeX const& gridx,
            IdxRangeVx const& gridv,
            RhsType const type,
            double extent,
            double stiffness,
            double amplitude,
            double density,
            double temperature);

    KrookSourceConstant(KrookSourceConstant&&) = default;

    ~KrookSourceConstant() override = default;

    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, double dt) const override;
};
```


