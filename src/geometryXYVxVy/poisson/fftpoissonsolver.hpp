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
            SplineXYBuilder const& spline_xy_builder,
            SplineXYEvaluator const& spline_xy_evaluator,
            SplineVxVyBuilder const& spline_vxvy_builder,
            SplineVxVyEvaluator const& spline_vxvy_evaluator);

    ~FftPoissonSolver() override = default;

    void operator()(
            DSpanXY electrostatic_potential,
            DSpanXY electric_field_x,
            DSpanXY electric_field_y,
            DViewSpXYVxVy allfdistribu) const override;
};
