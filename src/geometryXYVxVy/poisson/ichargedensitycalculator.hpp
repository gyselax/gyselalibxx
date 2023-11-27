// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

class IChargeDensityCalculator
{
public:
    IChargeDensityCalculator();

    void operator()(DSpanXY rho, DViewSpXYVxVy allfdistribu) const = 0;
};
