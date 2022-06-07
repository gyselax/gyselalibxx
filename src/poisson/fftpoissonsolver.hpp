// SPDX-License-Identifier: MIT

#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <fft.hpp>
#include <ifft.hpp>

#include "chargedensitycalculator.hpp"
#include "electricfield.hpp"
#include "ipoissonsolver.hpp"

class FftPoissonSolver : public IPoissonSolver
{
    IFourierTransform<RDimX> const& m_fft;

    IInverseFourierTransform<RDimX> const& m_ifft;

    SplineXBuilder const& m_spline_x_builder;

    SplineEvaluator<BSplinesX> const& m_spline_x_evaluator;

    ChargeDensityCalculator compute_rho;

    ElectricField m_electric_field;

public:
    FftPoissonSolver(
            IFourierTransform<RDimX> const& fft,
            IInverseFourierTransform<RDimX> const& ifft,
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator,
            SplineVxBuilder const& spline_vx_builder,
            SplineEvaluator<BSplinesVx> const& spline_vx_evaluator);

    ~FftPoissonSolver() override = default;

    void operator()(DSpanX electrostatic_potential, DSpanX electric_field, DViewSpXVx allfdistribu)
            const override;
};
