#pragma once

#include <fft.h>
#include <geometry.h>
#include <ifft.h>

#include "iefieldsolver.h"

class EfieldFftSolver : public IEfieldSolver
{
    IFourierTransform<Dim::X> const& m_fft;

    IInverseFourierTransform<Dim::X> const& m_ifft;

public:
    EfieldFftSolver(
            IFourierTransform<Dim::X> const& fft,
            IInverseFourierTransform<Dim::X> const& ifft);

    DSpanX operator()(DSpanX ex, DViewXVx fdistribu) const override;
};
