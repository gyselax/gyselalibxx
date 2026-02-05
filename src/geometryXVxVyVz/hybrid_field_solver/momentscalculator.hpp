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
    Quadrature<IdxRangeVxVyVz, IdxRangeXVxVyVz> m_quadrature;
    Quadrature<IdxRangeVz, IdxRangeXVxVyVz> m_quadrature_1D;

public:
    /**
     * @brief Create a ChargeDensityCalculator object.
     * @param[in] coeffs
     *            The coefficients of the quadrature.
     */
    explicit MomentsCalculator(DConstFieldVxVyVz coeffs, DConstFieldVz coeffs1D);

    /**
     * @brief Computes the charge density rho from the distribution function.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     */
    void operator()(DFieldX rho, DConstFieldSpVxVyVzX allfdistribu) const final;

    /**
     * @brief Computes the kinetic energy density rho from the distribution function.
     * @param[in, out] kinetic
     * @param[in] allfdistribu 
     */
    void operator()(DFieldX kinetic, DConstFieldSpVxVyVzX allfdistribu, char ki) const final;

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
    virtual void operator()(DFieldX mean_current_x, DFieldX mean_current_y, DFieldX mean_current_z, DFieldX rho, 
                                                 DConstFieldSpVxVyVzX allfdistribu) const final;

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
    virtual void operator()(DFieldX momentum_x, DFieldX momentum_y, DFieldX momentum_z, 
                                                 DConstFieldSpVxVyVzX allfdistribu) const final;

    /**
     * @brief Computes the charge density rho for each distribution function.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     */
    void operator()(DFieldSpX rho, DConstFieldSpVxVyVzX allfdistribu) const final;

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
    virtual void operator()(DFieldSpX mean_current_x, DFieldSpX mean_current_y, DFieldSpX mean_current_z, DFieldSpX rho, 
                                                 DConstFieldSpVxVyVzX allfdistribu) const final;

    /**
     * @brief Computes the density for each distribution function by intergrating for one velocity direction.
     * @param[in, out] rho
     * @param[in] allfdistribu 
     */
    void operator()(DFieldSpXVxVy rho, DConstFieldSpVxVyVzX allfdistribu) const final;

    
    void operator()(DFieldSpX parallel_temperature, DFieldSpX perpendicular_temperature, DConstFieldSpX rho, 
                    DConstFieldX Bx, DConstFieldX By, DConstFieldX Bz, 
                    DConstFieldSpX ux, DConstFieldSpX uy, DConstFieldSpX uz, 
                    DConstFieldSpVxVyVzX allfdistribu) const final;
};
