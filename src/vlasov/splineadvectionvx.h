#pragma once

#include <geometry.h>

#include "iadvectionvx.h"

class BoundaryValue;

class SplineAdvectionVx : public IAdvectionVx
{
private:
    const BSplinesVx& m_vx_spline_basis;

    const SplineVxBuilder& m_spline_vx_builder;

    const BoundaryValue& m_bc_left;

    const BoundaryValue& m_bc_right;

public:
    SplineAdvectionVx(const BSplinesVx& bspl, const SplineVxBuilder& spl_interp);

    SplineAdvectionVx(
            const BSplinesVx& bspl,
            const SplineVxBuilder& spl_interp,
            const BoundaryValue& bc_left,
            const BoundaryValue& bc_right);

    DSpanXVx operator()(DSpanXVx fdistribu, DViewX efield, double sqrt_me_on_mspecies, double dt)
            const override;
};
