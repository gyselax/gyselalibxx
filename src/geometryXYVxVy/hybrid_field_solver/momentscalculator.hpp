// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "geometry.hpp"
#include "imomentscalculator.hpp"
#include "quadrature.hpp"

/**
 * @brief A class which computes moments density with Kokkos.
 *
 * A class which computes charges density by solving the equation:
 * @f$ \int_{vx} \int_{vy} q_s v^r f_s(x,y,vx,vy) dvx dvy @f$
 * where @f$ q_s @f$ is the charge of the species @f$ s @f$ and
 * @f$ f_s(x,y,vx,vy) @f$ is the distribution function.
 */
class MomentsCalculator : public IMomentsCalculator
{
private:
    Quadrature<IdxRangeVxVy, IdxRangeXYVxVy> m_quadrature;

public:
    /**
     * @brief Create a ChargeDensityCalculator object.
     * @param[in] coeffs
     *            The coefficients of the quadrature.
     */
    explicit MomentsCalculator(DConstFieldVxVy coeffs);

    /**
     * @brief Computes the charge density rho from the distribution function.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     */
    void operator()(DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const final;

    /**
     * @brief Computes the kinetic energy density rho from the distribution function.
     * @param[in, out] kinetic
     * @param[in] allfdistribu 
     */
    void operator()(DFieldXY kinetic, DConstFieldSpVxVyXY allfdistribu, char ki) const final;

    /**
     * Calculate the mean velocity in direction x and y from the distribution function.
     *
     * Calculate the mean velocity by calculating the spline representation of slices
     * of the distribution function at each spatial point along the velocity direction.
     * This representation is multiplied with vx, vy and integrated, then multiplied by the charge,
     * and finally devided by the charge density.
     *
     * @param[out] mean_velocity_x The mean current in direction x.
     * @param[out] rho The charge density.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(DFieldXY mean_current_x, DFieldXY mean_current_y, DFieldXY rho, 
                                                 DConstFieldSpVxVyXY allfdistribu) const final;

    /**
     * Calculate the momentum in direction x and y from the distribution function.
     *
     * Calculate the momentum by calculating the spline representation of slices
     * of the distribution function at each spatial point along the velocity direction.
     * This representation is multiplied with vx, vy and integrated, and then multiplied by the mass,
     *
     * @param[out] mean_velocity_x The mean current in direction x.
     * @param[out] rho The charge density.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(DFieldXY momentum_x, DFieldXY momentum_y, 
                                                 DConstFieldSpVxVyXY allfdistribu) const final;

    /**
     * @brief Computes the charge density rho for each distribution function.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     */
    void operator()(DFieldSpXY rho, DConstFieldSpVxVyXY allfdistribu) const final;

    /**
     * Calculate the mean velocity in direction x and y for each distribution function.
     *
     * Calculate the mean velocity by calculating the spline representation of slices
     * of the distribution function at each spatial point along the velocity direction.
     * This representation is multiplied with vx, vy and integrated, then multiplied by the charge,
     * and finally divied by the charge density of the same species. 
     *
     * @param[out] mean_current_x The mean velocity in direction x.
     * @param[out] mean_current_y The mean velocity in direction y.
     * @param[out] rho The charge density.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(DFieldSpXY mean_current_x, DFieldSpXY mean_current_y, DFieldSpXY rho, 
                                                 DConstFieldSpVxVyXY allfdistribu) const final;
};
