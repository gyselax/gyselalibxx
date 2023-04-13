// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "iboltzmannsolver.hpp"

template <class Geometry, class DDimX>
class IAdvectionSpatial;
template <class Geometry, class DDimV>
class IAdvectionVelocity;

class SplitVlasovSolver : public IBoltzmannSolver
{
    IAdvectionSpatial<GeometryXVx, IDimX> const& m_advec_x;

    IAdvectionVelocity<GeometryXVx, IDimVx> const& m_advec_vx;

public:
    SplitVlasovSolver(
            IAdvectionSpatial<GeometryXVx, IDimX> const& advec_x,
            IAdvectionVelocity<GeometryXVx, IDimVx> const& m_advec_vx);

    ~SplitVlasovSolver() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_field, double dt) const override;
};
