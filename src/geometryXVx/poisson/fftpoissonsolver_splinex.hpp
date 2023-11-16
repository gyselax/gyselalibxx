// SPDX-License-Identifier: MIT

#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include "chargedensitycalculator_splineless.hpp"
#include "electricfield.hpp"
#include "ipoissonsolver.hpp"

class FftPoissonSolverSplineX : public IPoissonSolver
{
    ChargeDensityCalculatorSplineless m_compute_rho;
    ElectricField m_electric_field;

public:
    FftPoissonSolverSplineX(
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator,
            Quadrature<IDimVx> const& quad);

    ~FftPoissonSolverSplineX() override = default;

    void operator()(DSpanX electrostatic_potential, DSpanX electric_field, DViewSpXVx allfdistribu)
            const override;
};
