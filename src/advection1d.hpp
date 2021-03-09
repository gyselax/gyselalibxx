#pragma once
#include "null_boundary_value.hpp"
#include "spline_interpolator_1d.h"

// TODO: Move
/// A type for a 1D view of the memory
using View = std::experimental::basic_mdspan<double,
      std::experimental::extents<std::experimental::dynamic_extent>>;

class Advection1D
{
    public:
        Advection1D(const BSplines& bspl, const Spline_interpolator_1D& spl_interp,
                    double dt);
        Advection1D(const BSplines& bspl, const Spline_interpolator_1D& spl_interp,
                    double dt,
                    const BoundaryValue& bc_left,
                    const BoundaryValue& bc_right);
        void operator()(View& current_values, View&& velocity) const;
    private:
        const BSplines& m_bspl;
        const Spline_interpolator_1D& m_spline_interpolator;
        double dt;
        const BoundaryValue& m_bc_left;
        const BoundaryValue& m_bc_right;
};
