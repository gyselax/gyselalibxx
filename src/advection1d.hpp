#pragma once
#include "selalib_CXX/splines/spline_interpolator_1d.h"

// TODO: Move
/// A type for a 1D view of the memory
using View = std::experimental::basic_mdspan<double,
      std::experimental::extents<std::experimental::dynamic_extent>>;

class Advection1D
{
    public:
        Advection1D(const BSplines& bspl, const Spline_interpolator_1D& spl_interp,
                    double dt,
                    double (*bc_left)(double) = nullptr,
                    double (*bc_right)(double) = nullptr);
        void step(View& current_values, View&& velocity) const;
    private:
        const BSplines& m_bspl;
        const Spline_interpolator_1D& m_spline_interpolator;
        double (*m_bc_left)(double);
        double (*m_bc_right)(double);
        double dt;
};
