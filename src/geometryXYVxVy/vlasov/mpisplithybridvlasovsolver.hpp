// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "iadvectionvx.hpp"
#include "iadvectionrot2d.hpp"
#include "iadvectionx.hpp"
#include "ihybridvlasovsolver.hpp"
#include "mpitransposealltoall.hpp"

/**
 * @brief A class that solves a Vlasov equation using Strang's splitting
 * on an MPI distributed mesh.
 *
 * The Vlasov equation is split between advection equations in space and velocity.
 * The space advections involves the advections in X and Y directions. 
 * The velocity advections are done with the shifted exact splittings for rotations.
 */
class MpiSplitHybridVlasovSolver : public IHybridVlasovSolver
{
    /// Advection operator in the x direction
    IAdvectionSpatial<GeometryVxVyXY, GridX> const& m_advec_x;
    /// Advection operator in the y direction
    IAdvectionSpatial<GeometryVxVyXY, GridY> const& m_advec_y;
    /// Advection operator in the vx direction
    IAdvectionVelocity<GeometryXYVxVy, GridVx> const& m_advec_vx;
    /// Advection operator in the vy direction
    IAdvectionVelocity<GeometryXYVxVy, GridVy> const& m_advec_vy;
    /// Advection operator in the vx direction inside of the exact splitting
    IAdvectionVelocityRot2D<GeometryXYVxVy, GridVx> const& m_advec_2d_rot_vx;
    /// Advection operator in the vy direction inside of the exact splitting
    IAdvectionVelocityRot2D<GeometryXYVxVy, GridVy> const& m_advec_2d_rot_vy;
    /// MPI transpose operator
    MPITransposeAllToAll<X2DSplit, V2DSplit> const& m_transpose;

public:
    /**
     * @brief Creates an instance of the split vlasov solver class.
     * @param[in] advec_x An advection operator along the x direction.
     * @param[in] advec_y An advection operator along the y direction.
     * @param[in] advec_vx An advection operator along the vx direction.
     * @param[in] advec_vy An advection operator along the vy direction.
     * @param[in] advec_2d_rot_vx An advection operator along the vx direction inside of the exact splitting.
     * @param[in] advec_2d_rot_vy An advection operator along the vy direction inside of the exact splitting.
     * @param[in] transpose A MPI transpose operator to move between layouts.
     */
    MpiSplitHybridVlasovSolver(
            IAdvectionSpatial<GeometryVxVyXY, GridX> const& advec_x,
            IAdvectionSpatial<GeometryVxVyXY, GridY> const& advec_y,
            IAdvectionVelocity<GeometryXYVxVy, GridVx> const& advec_vx,
            IAdvectionVelocity<GeometryXYVxVy, GridVy> const& advec_vy,
            IAdvectionVelocityRot2D<GeometryXYVxVy, GridVx> const& advec_2d_rot_vx,
            IAdvectionVelocityRot2D<GeometryXYVxVy, GridVy> const& advec_2d_rot_vy,
            MPITransposeAllToAll<X2DSplit, V2DSplit> const& transpose);

    ~MpiSplitHybridVlasovSolver() override = default;

    /**
     * @brief Solves a Vlasov equation on a timestep dt.
     *
     * @param[in, out] allfdistribu On input : the initial value of the distribution function.
     *                              On output : the value of the distribution function after solving 
     *                              the Vlasov equation.
     * @param[in] magnetic_field_z The magnetic field in the z direction computed at all spatial positions. 
     * @param[in] frame_shift_x The velocity frame shift in x direction of the exact splitting.
     * @param[in] frame_shift_y The velocity frame shift in y direction of the exact splitting.
     * @param[in] dt The timestep. 
     *
     * @return The distribution function after solving the Vlasov equation.
     */
    DFieldSpVxVyXY operator()(
            DFieldSpVxVyXY allfdistribu,
            DFieldXY magnetic_field_z,
            DFieldXY frame_shift_x,
            DFieldXY frame_shift_y,
            double dt) const override;

    DFieldSpVxVyXY operator()(
            DFieldSpVxVyXY allfdistribu,
            double dt) const override;
};
