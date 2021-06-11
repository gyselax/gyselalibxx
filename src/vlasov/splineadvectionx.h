#pragma once

#include "blockview.h"
#include "boundary_value.h"
#include "iadvectionx.h"
#include "spline_builder_1d.h"

class SplineAdvectionX : public IAdvectionX
{
private:
    const deprecated::BSplines& m_x_spline_basis;

    const deprecated::SplineBuilder1D& m_spline_interpolator;

    const BoundaryValue& m_bc_left;

    const BoundaryValue& m_bc_right;

public:
    SplineAdvectionX(
            const deprecated::BSplines& bspl,
            const deprecated::SplineBuilder1D& spl_interp);

    SplineAdvectionX(
            const deprecated::BSplines& bspl,
            const deprecated::SplineBuilder1D& spl_interp,
            const BoundaryValue& bc_left,
            const BoundaryValue& bc_right);

    DBlockSpanXVx operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt) const override;
};
