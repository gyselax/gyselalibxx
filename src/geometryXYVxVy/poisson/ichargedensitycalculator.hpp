// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

/**
 * @brief A class which computes charges density.
 *
 * A class which computes charges density by solving the equation:
 * @f$ \int_{v} q_s f_s(x,v) dv @f$
 * where @f$ q_s @f$ is the charge of the species @f$ s @f$ and
 * @f$ f_s(x,v) @f$ is the distribution function.
 */
class IChargeDensityCalculator
{
public:
    /**
     * Calculate the charge density rho from the distribution function.
     *
     * Calculate the charge density by calculating the spline representation of slices
     * of the distribution function at each spatial point along the velocity direction.
     * This representation is then integrated and multiplied by the charge to find the
     * charge density.
     *
     * @param[out] rho The charge density.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(DSpanXY rho, DViewSpXYVxVy allfdistribu) const = 0;
};
