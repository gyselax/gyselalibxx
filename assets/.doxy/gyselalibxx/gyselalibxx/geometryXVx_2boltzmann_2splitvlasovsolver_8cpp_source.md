

# File splitvlasovsolver.cpp

[**File List**](files.md) **>** [**boltzmann**](dir_7559acab695a99e26dbd57f46ed1b0cd.md) **>** [**splitvlasovsolver.cpp**](geometryXVx_2boltzmann_2splitvlasovsolver_8cpp.md)

[Go to the documentation of this file](geometryXVx_2boltzmann_2splitvlasovsolver_8cpp.md)


```C++
// SPDX-License-Identifier: MIT

#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "splitvlasovsolver.hpp"

SplitVlasovSolver::SplitVlasovSolver(
        IAdvectionSpatial<GeometryXVx, GridX> const& advec_x,
        IAdvectionVelocity<GeometryXVx, GridVx> const& advec_vx)
    : m_advec_x(advec_x)
    , m_advec_vx(advec_vx)
{
}

DFieldSpXVx SplitVlasovSolver::operator()(
        DFieldSpXVx const allfdistribu,
        DConstFieldX const electric_field,
        double const dt) const
{
    m_advec_x(allfdistribu, dt / 2);
    m_advec_vx(allfdistribu, electric_field, dt);
    m_advec_x(allfdistribu, dt / 2);
    return allfdistribu;
}
```


