#pragma once

#include <sll/spline_evaluator.h>

#include <fft.h>
#include <geometry.h>
#include <ifft.h>

#include "iefieldsolver.h"

class EfieldFftSolver : public IEfieldSolver
{
    IFourierTransform<Dim::X> const& m_fft;

    IInverseFourierTransform<Dim::X> const& m_ifft;

    BSplinesX const& m_spline_x_basis;

    SplineXBuilder const& m_spline_x_builder;

    SplineEvaluator<BSplinesX> m_spline_x_evaluator;

    BSplinesVx const& m_spline_vx_basis;

    SplineVxBuilder const& m_spline_vx_builder;

    SplineEvaluator<BSplinesVx> m_spline_vx_evaluator;

    std::vector<double> m_derivs_vxmin_data;

    Span1D<double> m_derivs_vxmin;

    std::vector<double> m_derivs_vxmax_data;

    Span1D<double> m_derivs_vxmax;

public:
    EfieldFftSolver(
            IFourierTransform<Dim::X> const& fft,
            IInverseFourierTransform<Dim::X> const& ifft,
            BSplinesX const& bsplines_x,
            SplineXBuilder const& spline_x_builder,
            BSplinesVx const& bsplines_vx,
            SplineVxBuilder const& spline_vx_builder);

    DSpanX operator()(DSpanX ex, DViewXVx fdistribu) const override;
};
