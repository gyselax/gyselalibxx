// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "iadvectionvx.hpp"
#include "iadvectionrot3d.hpp"
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
    IAdvectionSpatial<GeometryVxVyVzX, GridX> const& m_advec_x;
    /// Advection operator in the vx direction
    IAdvectionVelocity<GeometryXVxVyVz, GridVx> const& m_advec_vx;
    /// Advection operator in the vy direction
    IAdvectionVelocity<GeometryXVxVyVz, GridVy> const& m_advec_vy;
    /// Advection operator in the vz direction
    IAdvectionVelocity<GeometryXVxVyVz, GridVz> const& m_advec_vz;
    /// Advection operators of the exact splitting
    IAdvectionVelocityRot3D<GeometryXVxVyVz, GridVx, GridVy, GridVz> const& m_advec_3d_rot;
    /// MPI transpose operator
    MPITransposeAllToAll<X1DSplit, V3DSplit> const& m_transpose;

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
            IAdvectionSpatial<GeometryVxVyVzX, GridX> const& advec_x,
            IAdvectionVelocity<GeometryXVxVyVz, GridVx> const& advec_vx,
            IAdvectionVelocity<GeometryXVxVyVz, GridVy> const& advec_vy,
            IAdvectionVelocity<GeometryXVxVyVz, GridVz> const& advec_vz,
            IAdvectionVelocityRot3D<GeometryXVxVyVz, GridVx, GridVy, GridVz> const& advec_3d_rot,
            MPITransposeAllToAll<X1DSplit, V3DSplit> const& transpose);

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
    DFieldSpVxVyVzX operator()(
            DFieldSpVxVyVzX allfdistribu,
            DFieldX frame_shift_x,
            DFieldX frame_shift_y,
            DFieldX frame_shift_z,
            DFieldX para_x,
            DFieldX para_y,
            DFieldX para_z,
            DFieldX B_x, DFieldX B_y, DFieldX B_z,
            double dt) const override;

    DFieldSpVxVyVzX operator()(
            DFieldSpVxVyVzX allfdistribu,
            double dt) const override;
};
