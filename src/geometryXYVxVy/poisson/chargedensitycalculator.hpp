// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

class ChargeDensityCalculator
{
    SplineVxVyBuilder const& m_spline_vxvy_builder;

    SplineVxVyEvaluator m_spline_vxvy_evaluator;

    int m_nbc_Vx;
    int m_nbc_Vy;
    int m_interp_dom_size_Vx;
    int m_interp_dom_size_Vy;

public:
    ChargeDensityCalculator(
            SplineVxVyBuilder const& spline_vxvy_builder,
            SplineVxVyEvaluator const& spline_vxvy_evaluator);

    void operator()(DSpanXY rho, DViewSpXYVxVy allfdistribu) const;
};
