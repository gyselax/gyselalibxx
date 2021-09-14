#pragma once

#include "geometry.h"
#include "iadvectionx.h"

class BoundaryValue;

class SplineAdvectionX : public IAdvectionX
{
private:
    const BSplinesX& m_x_spline_basis;

    const SplineXBuilder& m_spline_x_builder;

    const BoundaryValue& m_bc_left;

    const BoundaryValue& m_bc_right;

public:
    SplineAdvectionX(const BSplinesX& bspl, const SplineXBuilder& spl_interp);

    SplineAdvectionX(
            const BSplinesX& bspl,
            const SplineXBuilder& spl_interp,
            const BoundaryValue& bc_left,
            const BoundaryValue& bc_right);

    DSpanXVx operator()(DSpanXVx fdistribu, double mass_ratio, double dt) const override;
};
