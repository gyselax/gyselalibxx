

# File collisions\_inter.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**collisions\_inter.hpp**](collisions__inter_8hpp.md)

[Go to the documentation of this file](collisions__inter_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "irighthandside.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

class CollisionsInter : public IRightHandSide
{
private:
    double m_nustar0;
    DFieldMemSpX m_nustar_profile_alloc;
    DFieldSpX m_nustar_profile;

public:
    CollisionsInter(IdxRangeSpXVx const& mesh, double nustar0);

    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, double dt) const override;

    double get_nustar0() const;

    void get_derivative(DFieldSpXVx df, DConstFieldSpXVx allfdistribu) const;
};
```


