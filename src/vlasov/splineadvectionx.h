#pragma once

#include <sll/spline_evaluator.h>

#include <geometry.h>

#include "iadvectionx.h"

class BoundaryValue;
class SpeciesInformation;

class SplineAdvectionX : public IAdvectionX
{
private:
    const SplineXBuilder& m_spline_x_builder;

    SplineEvaluator<BSplinesX> m_spline_x_evaluator;

    SpeciesInformation const& m_species_info;

public:
    SplineAdvectionX(
            SpeciesInformation const& species,
            const BSplinesX& bspl,
            const SplineXBuilder& spl_interp);

    SplineAdvectionX(
            SpeciesInformation const& species,
            const BSplinesX& bspl,
            const SplineXBuilder& spl_interp,
            const BoundaryValue& bc_left,
            const BoundaryValue& bc_right);

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;
};
