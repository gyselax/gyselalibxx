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

public:
    EfieldFftSolver(
            IFourierTransform<Dim::X> const& fft,
            IInverseFourierTransform<Dim::X> const& ifft,
            BSplinesX const& bsplines_x,
            SplineXBuilder const& spline_x_builder);

    DSpanX operator()(DSpanX ex, DViewXVx fdistribu) const override;
};
