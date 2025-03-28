

# File splitvlasovsolver.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**vlasov**](dir_0a9688649b1824bbfb2c211b845ba732.md) **>** [**splitvlasovsolver.hpp**](geometryXYVxVy_2vlasov_2splitvlasovsolver_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2vlasov_2splitvlasovsolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "ivlasovsolver.hpp"

template <class Geometry, class GridX>
class IAdvectionSpatial;
template <class Geometry, class GridV>
class IAdvectionVelocity;

class SplitVlasovSolver : public IVlasovSolver
{
    IAdvectionSpatial<GeometryVxVyXY, GridX> const& m_advec_x;
    IAdvectionSpatial<GeometryVxVyXY, GridY> const& m_advec_y;

    IAdvectionVelocity<GeometryVxVyXY, GridVx> const& m_advec_vx;
    IAdvectionVelocity<GeometryVxVyXY, GridVy> const& m_advec_vy;

public:
    SplitVlasovSolver(
            IAdvectionSpatial<GeometryVxVyXY, GridX> const& advec_x,
            IAdvectionSpatial<GeometryVxVyXY, GridY> const& advec_y,
            IAdvectionVelocity<GeometryVxVyXY, GridVx> const& advec_vx,
            IAdvectionVelocity<GeometryVxVyXY, GridVy> const& advec_vy);

    ~SplitVlasovSolver() override = default;

    DFieldSpVxVyXY operator()(
            DFieldSpVxVyXY allfdistribu,
            DConstFieldXY electric_field_x,
            DConstFieldXY electric_field_y,
            double dt) const override;
};
```


