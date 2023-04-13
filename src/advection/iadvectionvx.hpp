// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

template <class Geometry, class DDimV>
class IAdvectionVelocity
{
public:
    virtual ~IAdvectionVelocity() = default;

    virtual ddc::ChunkSpan<double, typename Geometry::FdistribuDDom> operator()(
            ddc::ChunkSpan<double, typename Geometry::FdistribuDDom> allfdistribu,
            ddc::ChunkSpan<const double, typename Geometry::SpatialDDom> electrostatic_potential,
            double dt) const = 0;
};
