// SPDX-License-Identifier: MIT

#pragma once

#include "geometry_xyvxvy.hpp"
#include "ivlasovsolver.hpp"

template <class Geometry, class GridX>
class IAdvectionSpatial;
template <class Geometry, class GridV>
class IAdvectionVelocity;

/**
 * @brief A class that solves a Vlasov equation using Strang's splitting.
 *
 * The Vlasov equation is split between four advection equations 
 * along the X, Y, Vx and Vy directions. The splitting involves solving 
 * the advections in the X, Y, and Vx directions first on a time interval
 * of length dt/2, then the Vy-direction advection on a time dt, and
 * finally the X, Y, and Vx directions again in reverse order on dt/2.
 */
class SplitVlasovSolver : public IVlasovSolver
{
    /// Advection operator in the x direction
    IAdvectionSpatial<GeometryVxVyXY, GridX> const& m_advec_x;
    /// Advection operator in the y direction
    IAdvectionSpatial<GeometryVxVyXY, GridY> const& m_advec_y;

    /// Advection operator in the vx direction
    IAdvectionVelocity<GeometryVxVyXY, GridVx> const& m_advec_vx;
    /// Advection operator in the vy direction
    IAdvectionVelocity<GeometryVxVyXY, GridVy> const& m_advec_vy;

public:
    /**
     * @brief Creates an instance of the split vlasov solver class.
     * @param[in] advec_x An advection operator along the x direction.
     * @param[in] advec_y An advection operator along the y direction.
     * @param[in] advec_vx An advection operator along the vx direction.
     * @param[in] advec_vy An advection operator along the vy direction.
     */
    SplitVlasovSolver(
            IAdvectionSpatial<GeometryVxVyXY, GridX> const& advec_x,
            IAdvectionSpatial<GeometryVxVyXY, GridY> const& advec_y,
            IAdvectionVelocity<GeometryVxVyXY, GridVx> const& advec_vx,
            IAdvectionVelocity<GeometryVxVyXY, GridVy> const& advec_vy);

    ~SplitVlasovSolver() override = default;

    /**
     * @brief Solves a Vlasov equation on a timestep dt.
     *
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving 
     *                              the Vlasov equation.
     * @param[in] electric_field The electric field computed at all spatial positions.
     * @param[in] dt The timestep. 
     *
     * @return The distribution function after solving the Vlasov equation.
     */
    DFieldSpVxVyXY operator()(
            DFieldSpVxVyXY allfdistribu,
            DVectorConstFieldXY electric_field,
            double dt) const override;
};
