#pragma once

#include "blockview.h"
#include "boundary_value.h"
#include "iadvectionx.h"
#include "spline_interpolator_1d.h"

class AdvectionX : public IAdvectionX
{
private:
    const BSplines& m_x_spline_basis;

    const Spline_interpolator_1D& m_spline_interpolator;

    const BoundaryValue& m_bc_left;

    const BoundaryValue& m_bc_right;

public:
    AdvectionX(const BSplines& bspl, const Spline_interpolator_1D& spl_interp);

    AdvectionX(
            const BSplines& bspl,
            const Spline_interpolator_1D& spl_interp,
            const BoundaryValue& bc_left,
            const BoundaryValue& bc_right);

    DBlockViewXVx operator()(DBlockViewXVx fdistribu, double mass_ratio, double dt) const override;
};
