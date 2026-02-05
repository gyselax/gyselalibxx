// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"

/**
 * @brief An operator which solves the nonlinear equation about magnetic field, pressure, and 
 * mean velocity using a fast Fourier transform.
 *
 * An operator which solves the nonlinear equation about magnetic field, pressure, and 
 * mean velocity using a fast Fourier transform on a periodic index range.
 * This operator only works for equidistant points.
 */
class IHybridFieldSolver
{
public:
    virtual ~IHybridFieldSolver() = default;

    /**
     * The operator which solves the equation using the method described by the class. TO BE DELETED.
     *
     * @param[out] pressure_field The pressure, the result of the nonlinear solver.
     * @param[in] pressure_field_old The pressure in last time step, the input of the nonlinear solver.
     * @param[in] pressure_field_mid The pressure in middle time step, the intermidiate value of the nonlinear solver.
     * @param[in] pressure_field_previous The pressure in previous iteration, the intermidiate value of the nonlinear solver.
     * @param[out] magnetic_field_z The magnetic field, the result of the nonlinear solver.
     * @param[in] magnetic_field_z_old The magnetic field in last time step, the input of the nonlinear solver.
     * @param[in] magnetic_field_z_mid The magnetic field in middle time step, the intermidiate value of the nonlinear solver.
     * @param[in] mean_velocity_x_mid The intermidiate value of the nonlinear solver, the average of ux during a time step.
     * @param[in] mean_velocity_y_mid The intermidiate value of the nonlinear solver, the average of uy during a time step.
     * @param[out] rhs_1 the intermidiate value of the nonlinear solver, and the velocity frame shift in vx.
     * @param[out] rhs_2 the intermidiate value of the nonlinear solver, and the velocity frame shift in vy.
     * @param[in] rhs_3 the intermidiate value of the nonlinear solver.
     * @param[in] rhs_4 the intermidiate value of the nonlinear solver.
     * @param[in] mean_velocity_x The mean velocity in vx in the last time step.
     * @param[in] mean_velocity_y The mean velocity in vy in the last time step.
     * @param[in] mean_velocity_x_each The mean velocity in vx in the last time step for each species ions.
     * @param[in] mean_velocity_x_each The mean velocity in vx in the last time step for each species ions.
     * @param[in] rho_each The charge density in the last time step for each species ions.
     * @param[in] gradx_magnetic The derivative along x of the magnetic field, the intermidiate value of the nonlinear solver.
     * @param[in] grady_magnetic The derivative along y of the magnetic field, the intermidiate value of the nonlinear solver.
     * @param[in] gradx_pressure The derivative along x of the pressure field, the intermidiate value of the nonlinear solver.
     * @param[in] grady_pressure The derivative along y of the pressure field, the intermidiate value of the nonlinear solver.
     * @param[in] rho The charge density for all species ions.
     * @param[in] allfdistribu The distribution function.
     * @param[in] dt The time step.
     */

    virtual void operator()(
            DFieldX magnetic_field_y, 
            DFieldX magnetic_field_y_old,               
            DFieldX magnetic_field_y_mid, 
            DFieldX magnetic_field_y_previous,
            DFieldX magnetic_field_z, 
            DFieldX magnetic_field_z_old, 
            DFieldX magnetic_field_z_mid, 
            DFieldX magnetic_field_z_previous,
            DFieldSpX u_old_x, 
            DFieldSpX u_old_y, 
            DFieldSpX u_old_z, 
            DFieldX u_bar_x, 
            DFieldX u_bar_y, 
            DFieldX u_bar_z, 
            DFieldX rho,
            DFieldSpX rho_each,
            DFieldX magnetic_field_x,
            DFieldX gradx_rho,
            DFieldX gradx_magnetic_field_y_mid,
            DFieldX gradx_magnetic_field_z_mid,
            DFieldX rhs_1,
            DFieldX rhs_2,
            DFieldX rhs_3,
            DFieldX rhs_5,
            DFieldX rhs_6,
            DFieldX Mxx, DFieldX Mxy, DFieldX Mxz,
            DFieldX Myx, DFieldX Myy, DFieldX Myz,
            DFieldX Mzx, DFieldX Mzy, DFieldX Mzz,
            DFieldX weighted_u_x, DFieldX weighted_u_y, DFieldX weighted_u_z,
            DFieldX weighted_p_para_x, DFieldX weighted_p_para_y, DFieldX weighted_p_para_z,
            DFieldX qx, DFieldX qy, DFieldX qz,
            DFieldX p_parallel_x, DFieldX p_parallel_y, DFieldX p_parallel_z,
            DConstFieldSpVxVyVzX const allfdistribu,
            double const electron_temperature,
            double const dt) const = 0;


    virtual void operator()(
            DFieldSpX multi_para_tem, 
            DFieldSpX multi_perp_tem, 
            DFieldX single_para_tem, 
            DFieldX single_perp_tem) const = 0;

};
