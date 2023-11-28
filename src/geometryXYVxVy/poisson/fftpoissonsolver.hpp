// SPDX-License-Identifier: MIT

#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include "chargedensitycalculator.hpp"
#include "electricfield.hpp"
#include "ipoissonsolver.hpp"

/**
 * @brief An operator which solves the Poisson equation using a fast
 * Fourier transform.
 *
 * An operator which solves the Poisson equation:
 * @f$ - \frac{d^2 \phi}{dx^2} = \rho @f$
 * using a fast Fourier transform on a periodic domain.
 * This operator only works for equidistant points.
 *
 * The electric field, @f$ \frac{d \phi}{dx} @f$ is calculated using
 * a spline interpolation implemented in ElectricField.
 */
class FftPoissonSolver : public IPoissonSolver
{
    IChargeDensityCalculator const& m_compute_rho;

    ElectricField m_electric_field;

public:
    /**
     * Construct the FftPoissonSolver operator.
     *
     * @param spline_xy_builder A spline builder which calculates the coefficients of a spline representation.
     * @param spline_xy_evaluator A spline evaluator which provides the value of a spline representation from its coefficients.
     * @param compute_rho The operator which calculates the charge density, the right hand side of the equation.
     */
    FftPoissonSolver(
            SplineXYBuilder const& spline_xy_builder,
            SplineXYEvaluator const& spline_xy_evaluator,
            IChargeDensityCalculator const& compute_rho);

    ~FftPoissonSolver() override = default;

    /**
     * The operator which solves the equation using the method described by the class.
     *
     * @param[out] electrostatic_potential The electrostatic potential, the result of the poisson solver.
     * @param[out] electric_field_x The x-component of the electric field, the gradient of the electrostatic potential.
     * @param[out] electric_field_y The y-component of the electric field, the gradient of the electrostatic potential.
     * @param[in] allfdistribu The distribution function.
     */
    void operator()(
            DSpanXY electrostatic_potential,
            DSpanXY electric_field_x,
            DSpanXY electric_field_y,
            DViewSpXYVxVy allfdistribu) const override;
};
