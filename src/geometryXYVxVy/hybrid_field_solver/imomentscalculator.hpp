// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "geometry.hpp"

/**
 * @brief A class which computes velocity moments of the distribution functions.
 *
 * A class which computes velocity moments of the distribution:
 * @f$ \int_{v} q_s v^r f_s(x,v) dv @f$
 * where @f$ q_s @f$ is the charge of the species @f$ s @f$ and
 * @f$ f_s(x,v) @f$ is the distribution function.
 */
class IMomentsCalculator
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
    virtual void operator()(DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const = 0;

    /**
     * Calculate the kinetic energy density from the distribution function.
     *
     * Calculate the kinetic energy density by calculating the spline representation of slices
     * of the distribution function at each spatial point along the velocity direction.
     * This representation is then multiplied by 0.5 v v, integrated, and multiplied by the mass to find the
     * kinetic energy density.
     *
     * @param[out] kinetic The kinetic energy density.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(DFieldXY kinetic, DConstFieldSpVxVyXY allfdistribu, char ki) const = 0;

    /**
     * Calculate the mean velocity in direction x and y from the distribution function.
     *
     * Calculate the mean velocity by calculating the spline representation of slices
     * of the distribution function at each spatial point along the velocity direction.
     * This representation is multiplied with vx, vy and integrated, then multiplied by the charge,
     * and finally devided by rho.
     *
     * @param[out] mean_velocity_x The mean velocity in direction x.
     * @param[out] mean_velocity_y The mean velocity in direction y.
     * @param[out] rho The charge density.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(DFieldXY mean_current_x, DFieldXY mean_current_y, DFieldXY rho, 
        DConstFieldSpVxVyXY allfdistribu) const = 0;


    /**
     * Calculate the 1st order velocity moments in direction x and y from the distribution function.
     *
     * Calculate the 1st order velocity moments by calculating the spline representation of slices
     * of the distribution function at each spatial point along the velocity direction.
     * This representation is multiplied with vx, vy and integrated, then multiplied by the mass.
     *
     * @param[out] momentum_x The momentum in direction x.
     * @param[out] momentum_y The momentum in direction y.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(DFieldXY momentum_x, DFieldXY momentum_y, 
        DConstFieldSpVxVyXY allfdistribu) const = 0;
    
    /**
     * Calculate the charge density rho for each distribution function.
     *
     * Calculate the charge density by calculating the spline representation of slices
     * of the distribution function at each spatial point along the velocity direction.
     * This representation is then integrated and multiplied by the charge to find the
     * charge density.
     *
     * @param[out] rho The charge density.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(DFieldSpXY rho, DConstFieldSpVxVyXY allfdistribu) const = 0;

    /**
     * Calculate the mean velocity in direction x and y for each distribution function.
     *
     * Calculate the mean velocity by calculating the spline representation of slices
     * of the distribution function at each spatial point along the velocity direction.
     * This representation is multiplied with vx, vy and integrated, then multiplied by the charge,
     * and finally devided by the rho of the same species. 
     *
     * @param[out] mean_velocity_x The mean velocity in direction x.
     * @param[out] mean_velocity_y The mean velocity in direction y.
     * @param[in] rho The charge density.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(DFieldSpXY mean_current_x, DFieldSpXY mean_current_y, DFieldSpXY rho, 
        DConstFieldSpVxVyXY allfdistribu) const = 0;
};
