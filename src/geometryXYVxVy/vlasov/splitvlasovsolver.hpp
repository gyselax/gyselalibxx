// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "ivlasovsolver.hpp"

template <class Geometry, class DDimX>
class IAdvectionSpatial;
template <class Geometry, class DDimV>
class IAdvectionVelocity;

class SplitVlasovSolver : public IVlasovSolver
{
    IAdvectionSpatial<GeometryXYVxVy, IDimX> const& m_advec_x;
    IAdvectionSpatial<GeometryXYVxVy, IDimY> const& m_advec_y;

    IAdvectionVelocity<GeometryXYVxVy, IDimVx> const& m_advec_vx;
    IAdvectionVelocity<GeometryXYVxVy, IDimVy> const& m_advec_vy;

public:
    SplitVlasovSolver(
            IAdvectionSpatial<GeometryXYVxVy, IDimX> const& advec_x,
            IAdvectionSpatial<GeometryXYVxVy, IDimY> const& advec_y,
            IAdvectionVelocity<GeometryXYVxVy, IDimVx> const& advec_vx,
            IAdvectionVelocity<GeometryXYVxVy, IDimVy> const& advec_vy);

    ~SplitVlasovSolver() override = default;

    DSpanSpXYVxVy operator()(
            DSpanSpXYVxVy allfdistribu,
            DViewXY electric_field_x,
            DViewXY electric_field_y,
            double dt) const override;
};
