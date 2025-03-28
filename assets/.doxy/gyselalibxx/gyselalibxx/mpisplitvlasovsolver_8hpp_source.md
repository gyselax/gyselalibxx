

# File mpisplitvlasovsolver.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**vlasov**](dir_0a9688649b1824bbfb2c211b845ba732.md) **>** [**mpisplitvlasovsolver.hpp**](mpisplitvlasovsolver_8hpp.md)

[Go to the documentation of this file](mpisplitvlasovsolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "ivlasovsolver.hpp"
#include "mpitransposealltoall.hpp"

class MpiSplitVlasovSolver : public IVlasovSolver
{
    IAdvectionSpatial<GeometryVxVyXY, GridX> const& m_advec_x;
    IAdvectionSpatial<GeometryVxVyXY, GridY> const& m_advec_y;

    IAdvectionVelocity<GeometryXYVxVy, GridVx> const& m_advec_vx;
    IAdvectionVelocity<GeometryXYVxVy, GridVy> const& m_advec_vy;

    MPITransposeAllToAll<X2DSplit, V2DSplit> const& m_transpose;

public:
    MpiSplitVlasovSolver(
            IAdvectionSpatial<GeometryVxVyXY, GridX> const& advec_x,
            IAdvectionSpatial<GeometryVxVyXY, GridY> const& advec_y,
            IAdvectionVelocity<GeometryXYVxVy, GridVx> const& advec_vx,
            IAdvectionVelocity<GeometryXYVxVy, GridVy> const& advec_vy,
            MPITransposeAllToAll<X2DSplit, V2DSplit> const& transpose);

    ~MpiSplitVlasovSolver() override = default;

    DFieldSpVxVyXY operator()(
            DFieldSpVxVyXY allfdistribu,
            DConstFieldXY electric_field_x,
            DConstFieldXY electric_field_y,
            double dt) const override;
};
```


