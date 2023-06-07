// SPDX-License-Identifier: MIT

#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include "chargedensitycalculator.hpp"
#include "electricfield.hpp"
#include "ipoissonsolver.hpp"

class FftPoissonSolver : public IPoissonSolver
{
    ChargeDensityCalculator m_compute_rho;

    ElectricField m_electric_field;

public:
    FftPoissonSolver(
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator,
            SplineVxBuilder const& spline_vx_builder,
            SplineEvaluator<BSplinesVx> const& spline_vx_evaluator);

    ~FftPoissonSolver() override = default;

    void operator()(DSpanX electrostatic_potential, DSpanX electric_field, DViewSpXVx allfdistribu)
            const override;
};
