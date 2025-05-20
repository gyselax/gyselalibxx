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
     * @param[out] electrostatic_potential The electrostatic potential, the result of the poisson solver.
     * @param[out] electric_field_x The x-component of the electric field, the gradient of the electrostatic potential.
     * @param[out] electric_field_y The y-component of the electric field, the gradient of the electrostatic potential.
     * @param[in] allfdistribu The distribution function.
     */
    virtual void operator()(
            DFieldXY electrostatic_potential,
            DFieldXY electric_field_x,
            DFieldXY electric_field_y,
            DConstFieldSpVxVyXY allfdistribu) const = 0;

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
            DFieldXY pressure_field,
            DFieldXY pressure_field_old,
            DFieldXY pressure_field_mid,
            DFieldXY pressure_field_previous,
            DFieldXY magnetic_field_z,
            DFieldXY magnetic_field_z_old,
            DFieldXY magnetic_field_z_mid,
            DFieldXY magnetic_field_z_previous,
            DFieldXY mean_velocity_x_mid,
            DFieldXY mean_velocity_y_mid,
            DFieldXY rhs_1,
            DFieldXY rhs_2,
            DFieldXY rhs_3,
            DFieldXY rhs_4,
            DFieldXY mean_velocity_x,
            DFieldXY mean_velocity_y,
            DFieldSpXY mean_velocity_x_each,
            DFieldSpXY mean_velocity_y_each,
            DFieldSpXY rho_each,
            DFieldXY gradx_magnetic,
            DFieldXY grady_magnetic,
            DFieldXY gradx_pressure,
            DFieldXY grady_pressure,
            DFieldXY rho,
            DConstFieldSpVxVyXY const allfdistribu,
            double const dt) const = 0;
};
