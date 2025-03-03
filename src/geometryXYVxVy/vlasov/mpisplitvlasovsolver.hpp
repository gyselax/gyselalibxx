// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "ivlasovsolver.hpp"
#include "mpitransposealltoall.hpp"

/**
 * @brief A class that solves a Vlasov equation using Strang's splitting
 * on an MPI distributed mesh.
 *
 * The Vlasov equation is split between four advection equations 
 * along the X, Y, Vx and Vy directions. The splitting involves solving 
 * the advections in the X, Y, and Vx directions first on a time interval
 * of length dt/2, then the Vy-direction advection on a time dt, and
 * finally the X, Y, and Vx directions again in reverse order on dt/2.
 */
class MpiSplitVlasovSolver : public IVlasovSolver
{
    /// Advection operator in the x direction
    IAdvectionSpatial<GeometryVxVyXY, GridX> const& m_advec_x;
    /// Advection operator in the y direction
    IAdvectionSpatial<GeometryVxVyXY, GridY> const& m_advec_y;

    /// Advection operator in the vx direction
    IAdvectionVelocity<GeometryXYVxVy, GridVx> const& m_advec_vx;
    /// Advection operator in the vy direction
    IAdvectionVelocity<GeometryXYVxVy, GridVy> const& m_advec_vy;

    /// MPI transpose operator
    MPITransposeAllToAll<X2DSplit, V2DSplit> const& m_transpose;

public:
    /**
     * @brief Creates an instance of the split vlasov solver class.
     * @param[in] advec_x An advection operator along the x direction.
     * @param[in] advec_y An advection operator along the y direction.
     * @param[in] advec_vx An advection operator along the vx direction.
     * @param[in] advec_vy An advection operator along the vy direction.
     */
    MpiSplitVlasovSolver(
            IAdvectionSpatial<GeometryVxVyXY, GridX> const& advec_x,
            IAdvectionSpatial<GeometryVxVyXY, GridY> const& advec_y,
            IAdvectionVelocity<GeometryXYVxVy, GridVx> const& advec_vx,
            IAdvectionVelocity<GeometryXYVxVy, GridVy> const& advec_vy,
            MPITransposeAllToAll<X2DSplit, V2DSplit> const& transpose);

    ~MpiSplitVlasovSolver() override = default;

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
    DFieldSpVxVyXY operator()(
            DFieldSpVxVyXY allfdistribu,
            DConstFieldXY electric_field_x,
            DConstFieldXY electric_field_y,
            double dt) const override;
};
