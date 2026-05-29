

# File splitvlasovsolver.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**vlasov**](dir_0a9688649b1824bbfb2c211b845ba732.md) **>** [**splitvlasovsolver.hpp**](geometryXYVxVy_2vlasov_2splitvlasovsolver_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2vlasov_2splitvlasovsolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry_xyvxvy.hpp"
#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "ivlasovsolver.hpp"

class SplitVlasovSolver : public IVlasovSolver
{
    IAdvectionSpatial<GeometryVxVyXY, GridX, Real> const& m_advec_x;
    IAdvectionSpatial<GeometryVxVyXY, GridY, Real> const& m_advec_y;

    IAdvectionVelocity<GeometryVxVyXY, GridVx, Real> const& m_advec_vx;
    IAdvectionVelocity<GeometryVxVyXY, GridVy, Real> const& m_advec_vy;

public:
    SplitVlasovSolver(
            IAdvectionSpatial<GeometryVxVyXY, GridX, Real> const& advec_x,
            IAdvectionSpatial<GeometryVxVyXY, GridY, Real> const& advec_y,
            IAdvectionVelocity<GeometryVxVyXY, GridVx, Real> const& advec_vx,
            IAdvectionVelocity<GeometryVxVyXY, GridVy, Real> const& advec_vy);

    ~SplitVlasovSolver() override = default;

    DFieldSpVxVyXY operator()(
            DFieldSpVxVyXY allfdistribu,
            DVectorConstFieldXY electric_field,
            Real dt) const override;
};
```


