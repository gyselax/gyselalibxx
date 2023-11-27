// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

class IChargeDensityCalculator
{
public:
    IChargeDensityCalculator();

    void operator()(DSpanX rho, DViewSpXVx allfdistribu) const = 0;
};
