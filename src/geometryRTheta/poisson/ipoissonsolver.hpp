// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

/**
 * @brief Base class for Poisson solver.
 */
class IPoissonSolver
{
public:
    virtual ~IPoissonSolver() = default;

    /**
     * @brief Compute the electrical potential and
     * the electric field from the Poisson equation.
     *
     * @param[out] electrostatic_potential
     *      The solution of the Poisson equation.
     * @param[out] electric_field
     *      The electric field @f$E = -\nabla \phi@f$.
     * @param[in] allfdistribu
     *      The rhs of the Poisson equation.
     */
    virtual void operator()(
            DSpanRP electrostatic_potential,
            VectorFieldSpan<double, IDomainRP, NDTag<RDimX, RDimY>> electric_field,
            DViewRP allfdistribu) const = 0;
};
