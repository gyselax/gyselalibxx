// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

/**
 * @brief An abstract class representing a source for both the Boltzmann and neutral fluid equations.
**/
class IKineticFluidCoupling
{
public:
    virtual ~IKineticFluidCoupling() = default;

    /**
     * @brief Operator for applying the source term on the distribution function and neutral density.
     * 
     * The source @f$S@f$ acts on the distribution function following the 
     * evolution equation @f$df/dt = S_{n,N}(x) * S_v(x,v)@f$ where
     * @f$S_{n,N}(x) = n_N(x) n_e(x) K_i(x) - n_i(x) n_e(x) K_r(x)@f$
     * @f$S_v(x,v)@f$ is the sum of order 0 to 2 Hermite polynomials times a Maxwellian velocity distribution function.
     * 
     * The source @f$S@f$ acts on the neutral density following the
     * evolution equation @f$dn_N/dt = - S_n,N(x)@f$
     * 
     * @param[in, out] allfdistribu On input: the initial value of the distribution function.
     *                              On output: the value of the distribution function after solving 
     *                              the source evolution equation on one timestep.
     * @param[in, out] neutrals On input: the initial value of the neutral density.
     *                          On output: the value of the neutral density after solving
     *                          the source evolution equation on one timestep.
     * @param[in] dt The timestep.
     * 
     */
    virtual void operator()(DSpanSpXVx const allfdistribu, DSpanSpMX neutrals, double const dt)
            const = 0;
};
