#pragma once

#include <vector>

#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

#include "iadvectionvx.hpp"

class BoundaryValue;
class SpeciesInformation;

class SplineAdvectionVx : public IAdvectionVx
{
private:
    SplineXBuilder const& m_spline_x_builder;

    SplineEvaluator<BSplinesX> const& m_spline_x_evaluator;

    SplineVxBuilder const& m_spline_vx_builder;

    SplineEvaluator<BSplinesVx> const& m_spline_vx_evaluator;

    std::vector<double> m_derivs_vxmin_data;

    Span1D<double> m_derivs_vxmin;

    std::vector<double> m_derivs_vxmax_data;

    Span1D<double> m_derivs_vxmax;

    SpeciesInformation const& m_species_info;

public:
    SplineAdvectionVx(
            SpeciesInformation const& species_info,
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator,
            SplineVxBuilder const& spline_vx_builder,
            SplineEvaluator<BSplinesVx> const& spline_vx_evaluator);

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_potential, double dt)
            const override;
};
