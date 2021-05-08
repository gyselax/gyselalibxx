#pragma once

#include "block.h"
#include "boundary_value.h"
#include "spline_interpolator_1d.h"
#include "view.h"

class Advection1D
{
public:
    Advection1D(const BSplines& bspl, const Spline_interpolator_1D& spl_interp);

    Advection1D(
            const BSplines& bspl,
            const Spline_interpolator_1D& spl_interp,
            const BoundaryValue& bc_left,
            const BoundaryValue& bc_right);

    void operator()(DBlockView2D& fdistribu, double mass_ratio, double dt) const;

private:
    const BSplines& m_bspl;

    const Spline_interpolator_1D& m_spline_interpolator;

    const BoundaryValue& m_bc_left;

    const BoundaryValue& m_bc_right;
};
