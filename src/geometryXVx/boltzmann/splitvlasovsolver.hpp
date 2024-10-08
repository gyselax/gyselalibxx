// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "iboltzmannsolver.hpp"

/**A generic class for a spatial advection*/
template <class Geometry, class GridX>
class IAdvectionSpatial;
/**A generic class for a velocity advection*/
template <class Geometry, class GridV>
class IAdvectionVelocity;

/**
 * @brief A class that solves a Vlasov equation using Strang's splitting.
 *
 * The Vlasov equation is split between two advection equations 
 * along the X and Vx directions. The splitting involves solving 
 * the x-direction advection first on a time interval of length dt/2, 
 * then the vx-direction advection on a time dt, and then x-direction
 * again on dt/2.
 */
class SplitVlasovSolver : public IBoltzmannSolver
{
    /** Member advection operator in the x direction*/
    IAdvectionSpatial<GeometryXVx, GridX> const& m_advec_x;

    /** Member advection operator in the vx direction*/
    IAdvectionVelocity<GeometryXVx, GridVx> const& m_advec_vx;

public:
    /**
     * @brief Creates an instance of the split vlasov solver class.
     * @param[in] advec_x An advection operator along the x direction.
     * @param[in] advec_vx An advection operator along the vx direction.
     */
    SplitVlasovSolver(
            IAdvectionSpatial<GeometryXVx, GridX> const& advec_x,
            IAdvectionVelocity<GeometryXVx, GridVx> const& advec_vx);

    ~SplitVlasovSolver() override = default;

    /**
     * @brief Solves a Vlasov equation on a timestep dt.
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving 
     *                              the Vlasov equation.
     * @param[in] electric_field The electric field computed at all spatial positions. 
     * @param[in] dt The timestep. 
     * @return The distribution function after solving the Vlasov equation.
     */
    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, DConstFieldX electric_field, double dt)
            const override;
};
