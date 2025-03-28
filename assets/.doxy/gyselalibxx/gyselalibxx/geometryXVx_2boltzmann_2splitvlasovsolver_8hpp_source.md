

# File splitvlasovsolver.hpp

[**File List**](files.md) **>** [**boltzmann**](dir_7559acab695a99e26dbd57f46ed1b0cd.md) **>** [**splitvlasovsolver.hpp**](geometryXVx_2boltzmann_2splitvlasovsolver_8hpp.md)

[Go to the documentation of this file](geometryXVx_2boltzmann_2splitvlasovsolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "iboltzmannsolver.hpp"

template <class Geometry, class GridX>
class IAdvectionSpatial;
template <class Geometry, class GridV>
class IAdvectionVelocity;

class SplitVlasovSolver : public IBoltzmannSolver
{
    IAdvectionSpatial<GeometryXVx, GridX> const& m_advec_x;

    IAdvectionVelocity<GeometryXVx, GridVx> const& m_advec_vx;

public:
    SplitVlasovSolver(
            IAdvectionSpatial<GeometryXVx, GridX> const& advec_x,
            IAdvectionVelocity<GeometryXVx, GridVx> const& advec_vx);

    ~SplitVlasovSolver() override = default;

    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, DConstFieldX electric_field, double dt)
            const override;
};
```


