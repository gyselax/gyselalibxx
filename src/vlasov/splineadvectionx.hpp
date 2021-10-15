#pragma once

#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

#include "iadvectionx.hpp"

class BoundaryValue;
class SpeciesInformation;

class SplineAdvectionX : public IAdvectionX
{
private:
    const SplineXBuilder& m_spline_x_builder;

    SplineEvaluator<BSplinesX> const& m_spline_x_evaluator;

    SpeciesInformation const& m_species_info;

public:
    SplineAdvectionX(
            SpeciesInformation const& species_info,
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator);

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;
};
