// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

/**
 * @brief An enum class representing the type of a source.
 */
enum class RhsType { Source, Sink };

/**
 * @brief An abstract class representing a source in Boltzmann equation.
 */
class IRightHandSide
{
public:
    virtual ~IRightHandSide() = default;

    /**
     * @brief Operator for applying the source term on the distribution function.
     * 
     * The source @f$S@f$ acts on the distribution function following the 
     * evolution equation df/dt = S
     * @param[in, out] allfdistribu On input: the initial value of the distribution function.
     *                              On output: the value of the distribution function after solving 
     *                              the source evolution equation on one timestep.
     * @param[in] dt The timestep.
     * @return The distribution function after solving the source evolution equation.
     */
    virtual device_t<DSpanSpXVx> operator()(device_t<DSpanSpXVx> allfdistribu, double dt) const = 0;
};
