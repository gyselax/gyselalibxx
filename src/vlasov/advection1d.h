#pragma once

#include "block.h"
#include "boundary_value.h"
#include "iadvectionx.h"
#include "spline_interpolator_1d.h"

class Advection1D : public IAdvectionX
{
public:
    Advection1D(const BSplines& bspl, const Spline_interpolator_1D& spl_interp);

    Advection1D(
            const BSplines& bspl,
            const Spline_interpolator_1D& spl_interp,
            const BoundaryValue& bc_left,
            const BoundaryValue& bc_right);

    DBlockViewXVx& operator()(DBlockViewXVx& fdistribu, double mass_ratio, double dt)
            const override;

private:
    const BSplines& m_x_spline_basis;

    const Spline_interpolator_1D& m_spline_interpolator;

    const BoundaryValue& m_bc_left;

    const BoundaryValue& m_bc_right;
};
