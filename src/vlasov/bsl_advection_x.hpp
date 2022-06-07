// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "iadvectionx.hpp"

class IPreallocatableInterpolatorX;

class BslAdvectionX : public IAdvectionX
{
private:
    IPreallocatableInterpolatorX const& m_interpolator;

public:
    BslAdvectionX(IPreallocatableInterpolatorX const& interpolator);

    ~BslAdvectionX() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;
};
