// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for solving the hybrid plasma model with kinetic ions 
 *        and massless electrons.
 */
class IHybridSplitting
{
public:
    virtual ~IHybridSplitting() = default;

    /**
     * @brief Solves the hybrid model with kinetic ions and massless electrons.
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving
     *                              the hybrid model a given number of iterations.
     * @param[in] mean_velocity_x_each The mean velocity in vx for each species.
     * @param[in] mean_velocity_y_each The mean velocity in vy for each species.
     * @param[in] mean_velocity_x The mean velocity in vx for all species.
     * @param[in] momentum_x The momentum in vx for all species.
     * @param[in] momentum_y The momentum in vy for all species.
     * @param[in, out] magnetic_field_z The magnetic field.
     * @param[in, out] pressure The pressure field.
     * @param[in] rho_each The charge density for each species.
     * @param[in] rho The charge density for all species.
     * @param[out] kinetic The kinetic density for all species.
     * @param[in] dt The timestep.
     * @param[in] steps The number of iterations to be performed by the predictor-corrector.
     * @return The distribution function after solving the system.
     */
    virtual DFieldSpVxVyXY operator()(DFieldSpVxVyXY allfdistribu, 
                                      DFieldSpXY mean_velocity_x_each, 
                                      DFieldSpXY mean_velocity_y_each,
                                      DFieldXY mean_velocity_x,
                                      DFieldXY mean_velocity_y,
                                      DFieldXY momentum_x,
                                      DFieldXY momentum_y,
                                      DFieldXY magnetic_field_z,
                                      DFieldXY pressure,
                                      DFieldSpXY rho_each, 
                                      DFieldXY rho, 
                                      DFieldXY kinetic, 
                                      double dt, int steps = 1) const = 0;

};
