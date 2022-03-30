// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

class SpeciesInformation;

class ChargeDensityCalculator
{
    SpeciesInformation const& m_species_info;

    SplineVxBuilder const& m_spline_vx_builder;

    SplineEvaluator<BSplinesVx> m_spline_vx_evaluator;

    std::vector<double> m_derivs_vxmin_data;

    Span1D<double> m_derivs_vxmin;

    std::vector<double> m_derivs_vxmax_data;

    Span1D<double> m_derivs_vxmax;

public:
    ChargeDensityCalculator(
            SpeciesInformation const& species_info,
            SplineVxBuilder const& spline_vx_builder,
            SplineEvaluator<BSplinesVx> const& spline_vx_evaluator);

    void operator()(DSpanX const rho, DViewSpXVx const allfdistribu) const;
};
