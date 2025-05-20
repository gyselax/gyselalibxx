// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for solving the Vlasov equation in the hybrid model 
 * with kinetic ions and massless electrons by advections in space and velocity.
 */
class IHybridVlasovSolver
{
public:
    virtual ~IHybridVlasovSolver() = default;

    /**
     * @brief Solves a Vlasov equation on a timestep dt by velocity advections.
     *
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving 
     *                              the Vlasov equation.
     * @param[in] magnetic_field_z The magnetic field in the z direction computed at all spatial positions. 
     * @param[in] frame_shift_x The shift of the velocity frame in x direction.
     * @param[in] frame_shift_y The shift of the velocity frame in y direction. 
     * @param[in] dt The timestep. 
     *
     * @return The distribution function after solving the Vlasov equation.
     */
    virtual DFieldSpVxVyXY operator()(
            DFieldSpVxVyXY allfdistribu,
            DFieldXY magnetic_field_z,
            DFieldXY frame_shift_x,
            DFieldXY frame_shift_y,
            double dt) const = 0;

    /**
     * @brief Solves a Vlasov equation on a timestep dt by space advections.
     *
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving 
     *                              the Vlasov equation.
     * @param[in] dt The timestep. 
     *
     * @return The distribution function after solving the Vlasov equation.
     */
    virtual DFieldSpVxVyXY operator()(
            DFieldSpVxVyXY allfdistribu,
            double dt) const = 0;
};
