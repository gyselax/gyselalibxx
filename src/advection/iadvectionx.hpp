// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

template <class Geometry, class DDimX>
class IAdvectionSpatial
{
public:
    virtual ~IAdvectionSpatial() = default;

    virtual ddc::ChunkSpan<double, typename Geometry::FdistribuDDom> operator()(
            ddc::ChunkSpan<double, typename Geometry::FdistribuDDom> allfdistribu,
            double dt) const = 0;
};
