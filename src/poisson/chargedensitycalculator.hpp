// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

class ChargeDensityCalculator
{
    SplineVxBuilder const& m_spline_vx_builder;

    SplineEvaluator<BSplinesVx> m_spline_vx_evaluator;

    std::vector<double> m_derivs_vxmin_data;

    Span1D<double> m_derivs_vxmin;

    std::vector<double> m_derivs_vxmax_data;

    Span1D<double> m_derivs_vxmax;

public:
    ChargeDensityCalculator(
            SplineVxBuilder const& spline_vx_builder,
            SplineEvaluator<BSplinesVx> const& spline_vx_evaluator);

    void operator()(DSpanX rho, DViewSpXVx allfdistribu) const;
};
