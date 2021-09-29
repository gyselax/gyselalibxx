#pragma once

#include <sll/spline_evaluator.h>

#include <geometry.h>

#include "iadvectionx.h"

class BoundaryValue;

class SplineAdvectionX : public IAdvectionX
{
private:
    const BSplinesX& m_x_spline_basis;

    const SplineXBuilder& m_spline_x_builder;

    SplineEvaluator<BSplinesX> m_spline_x_evaluator;

public:
    SplineAdvectionX(const BSplinesX& bspl, const SplineXBuilder& spl_interp);

    SplineAdvectionX(
            const BSplinesX& bspl,
            const SplineXBuilder& spl_interp,
            const BoundaryValue& bc_left,
            const BoundaryValue& bc_right);

    DSpanXVx operator()(DSpanXVx fdistribu, double sqrt_me_on_mspecies, double dt) const override;
};
