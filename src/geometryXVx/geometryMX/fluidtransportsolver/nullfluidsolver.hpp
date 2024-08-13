// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "ifluidtransportsolver.hpp"

/**
 * @brief A dommy class that solves a fluid model.
 * The fluid model leaves the moments of the fluid species unchanged.
 */
class NullFluidSolver : public IFluidTransportSolver
{
public:
    /**
     * @brief The constructor for the class.
     *
     * @param[in] dom_fluidsp The moments index range on which the fluid species is defined.
     */
    NullFluidSolver(IdxRangeSp const& dom_fluidsp);

    ~NullFluidSolver() override = default;

    /**
     * @brief Solves a dummy fluid model on a timestep dt.
     * @param[in, out] fluid_moments On input : a field referencing the moments of the fluid species.
     *                               On output : a field referencing the moments of the fluid species 
     *                               updated after solving the dummy fluid model.
     * @param[in] allfdistribu A constant view referencing the distribution function.
     * @param[in] efield A constant view referencing the electric field.
     * @param[in] dt The timestep. 
     * @return a field referencing the fluid species after solving the dummy fluid model on one timestep.
     */
    DFieldSpMomX operator()(
            DFieldSpMomX fluid_moments,
            DConstFieldSpXVx allfdistribu,
            DConstFieldX efield,
            double dt) const override;
};
