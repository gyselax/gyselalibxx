

# File krook\_source\_adaptive.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**krook\_source\_adaptive.hpp**](krook__source__adaptive_8hpp.md)

[Go to the documentation of this file](krook__source__adaptive_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "geometry.hpp"
#include "irighthandside.hpp"

class KrookSourceAdaptive : public IRightHandSide
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
    KrookSourceAdaptive(
            IdxRangeX const& gridx,
            IdxRangeVx const& gridvx,
            RhsType const type,
            double extent,
            double stiffness,
            double amplitude,
            double density,
            double temperature);

    KrookSourceAdaptive(KrookSourceAdaptive&&) = default;

    ~KrookSourceAdaptive() override = default;

    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, double dt) const override;

public:
    void get_amplitudes(DFieldSpX amplitudes, DConstFieldSpXVx allfdistribu) const;

    void get_derivative(DFieldSpXVx df, DConstFieldSpXVx f, DConstFieldSpXVx f0) const;
};
```


