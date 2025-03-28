

# File splitvlasovsolver.cpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**vlasov**](dir_0a9688649b1824bbfb2c211b845ba732.md) **>** [**splitvlasovsolver.cpp**](geometryXYVxVy_2vlasov_2splitvlasovsolver_8cpp.md)

[Go to the documentation of this file](geometryXYVxVy_2vlasov_2splitvlasovsolver_8cpp.md)


```C++
// SPDX-License-Identifier: MIT

#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "splitvlasovsolver.hpp"

SplitVlasovSolver::SplitVlasovSolver(
        IAdvectionSpatial<GeometryVxVyXY, GridX> const& advec_x,
        IAdvectionSpatial<GeometryVxVyXY, GridY> const& advec_y,
        IAdvectionVelocity<GeometryVxVyXY, GridVx> const& advec_vx,
        IAdvectionVelocity<GeometryVxVyXY, GridVy> const& advec_vy)
    : m_advec_x(advec_x)
    , m_advec_y(advec_y)
    , m_advec_vx(advec_vx)
    , m_advec_vy(advec_vy)
{
}

DFieldSpVxVyXY SplitVlasovSolver::operator()(
        DFieldSpVxVyXY const allfdistribu,
        DConstFieldXY const electric_field_x,
        DConstFieldXY const electric_field_y,
        double const dt) const
{
    m_advec_x(allfdistribu, dt / 2);
    m_advec_y(allfdistribu, dt / 2);
    m_advec_vx(allfdistribu, electric_field_x, dt / 2);
    m_advec_vy(allfdistribu, electric_field_y, dt);
    m_advec_vx(allfdistribu, electric_field_x, dt / 2);
    m_advec_y(allfdistribu, dt / 2);
    m_advec_x(allfdistribu, dt / 2);

    return allfdistribu;
}
```


