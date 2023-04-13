// SPDX-License-Identifier: MIT

#pragma once

#include "iadvectionvx.hpp"

template <class FdistribuDDom, class SpatialDDom>
class NullAdvectionVelocity : public IAdvectionV<FdistribuDDom, SpatialDDom>
{
public:
    NullAdvectionVelocity() = default;

    ~NullAdvectionVelocity() override = default;

    ddc::ChunkSpan<double, FdistribuDDom> operator()(
            ddc::ChunkSpan<double, FdistribuDDom> allfdistribu,
            [[maybe_unused]] ddc::ChunkSpan<const double, SpatialDDom> electrostatic_potential,
            [[maybe_unused]] double dt) const override
    {
        return allfdistribu;
    }
};
