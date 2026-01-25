// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
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
    IAdvectionSpatial<GeometryVxVyVzX, GridX> const& m_advec_x;

    /// Advection operator in the vx direction
    IAdvectionVelocity<GeometryVxVyVzX, GridVx> const& m_advec_vx;
    /// Advection operator in the vy direction
    IAdvectionVelocity<GeometryVxVyVzX, GridVy> const& m_advec_vy;
    /// Advection operator in the vz direction
    IAdvectionVelocity<GeometryVxVyVzX, GridVz> const& m_advec_vz;

public:
    /**
     * @brief Creates an instance of the split vlasov solver class.
     * @param[in] advec_x An advection operator along the x direction.
     * @param[in] advec_y An advection operator along the y direction.
     * @param[in] advec_vx An advection operator along the vx direction.
     * @param[in] advec_vy An advection operator along the vy direction.
     */
    SplitVlasovSolver(
            IAdvectionSpatial<GeometryVxVyVzX, GridX> const& advec_x,
            IAdvectionVelocity<GeometryVxVyVzX, GridVx> const& advec_vx,
            IAdvectionVelocity<GeometryVxVyVzX, GridVy> const& advec_vy,
            IAdvectionVelocity<GeometryVxVyVzX, GridVz> const& advec_vz);

    ~SplitVlasovSolver() override = default;

    /**
     * @brief Solves a Vlasov equation on a timestep dt.
     *
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving 
     *                              the Vlasov equation.
     * @param[in] electric_field_x The electric field in the x direction computed at all spatial positions. 
     * @param[in] electric_field_y The electric field in the y direction computed at all spatial positions. 
     * @param[in] dt The timestep. 
     *
     * @return The distribution function after solving the Vlasov equation.
     */
    DFieldSpVxVyVzX operator()(
            DFieldSpVxVyVzX allfdistribu,
            DConstFieldX electric_field_x,
            DConstFieldX electric_field_y,
            DConstFieldX electric_field_z,
            double dt) const override;
};
