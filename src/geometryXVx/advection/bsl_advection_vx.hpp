// SPDX-License-Identifier: MIT

#pragma once

#include <vector>

#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

#include "iadvectionvx.hpp"

class IPreallocatableInterpolatorX;
class IPreallocatableInterpolatorVx;
class BoundaryValue;

class BslAdvectionVx : public IAdvectionVx
{
private:
    IPreallocatableInterpolatorVx const& m_interpolator_vx;

public:
    explicit BslAdvectionVx(IPreallocatableInterpolatorVx const& interpolator_vx);

    ~BslAdvectionVx() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_field, double dt) const override;
};
