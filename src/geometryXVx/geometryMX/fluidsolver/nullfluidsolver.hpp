// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "ifluidsolver.hpp"

/**
 * @brief A dommy class that solves a fluid model.
 * The fluid model leaves the moments of the fluid species unchanged.
 */
class NullFluidSolver : public IFluidSolver
{
public:
    /**
     * @brief The constructor for the class.
     *
     * @param[in] dom_fluidsp The moments domain on which the fluid species is defined.
     */
    NullFluidSolver(IDomainSp const& dom_fluidsp);

    ~NullFluidSolver() override = default;

    /**
     * @brief Solves a dummy fluid model on a timestep dt.
     * @param[in, out] fluid_moments On input : a span referencing the moments of the fluid species.
     *                               On output : a span referencing the moments of the fluid species 
     *                               updated after solving the dummy fluid model.
     * @param[in] allfdistribu A constant view referencing the distribution function.
     * @param[in] efield A constant view referencing the electric field.
     * @param[in] dt The timestep. 
     * @return a span referencing the fluid species after solving the dummy fluid model on one timestep.
     */
    DSpanSpMX operator()(DSpanSpMX fluid_moments, DViewSpXVx allfdistribu, DViewX efield, double dt)
            const override;
};