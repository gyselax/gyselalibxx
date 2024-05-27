// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

/**
 * A class which calculates the charge density.
 *
 * A class which calculates the charge density. This is then used as
 * the right hand side of the Quasi-Neutrality equation.
 */
class IChargeDensityCalculator
{
public:
    /**
     * Calculate the charge density rho from the distribution function.
     *
     * @param[out] rho The charge density.
     * @param[in] allfdistribu The distribution function.
     *
     * @return rho The charge density.
     */
    virtual DSpanX operator()(DSpanX rho, DViewSpXVx allfdistribu) const = 0;
};
