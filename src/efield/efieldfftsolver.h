#pragma once

#include "fftw.h"
#include "iefieldsolver.h"
#include "ifftw.h"

class EfieldFftSolver : public IEfieldSolver
{
    IFourierTransform<Dim::X> const& m_fft;

    IInverseFourierTransform<Dim::X> const& m_ifft;

public:
    EfieldFftSolver(
            IFourierTransform<Dim::X> const& fft,
            IInverseFourierTransform<Dim::X> const& ifft);

    DBlockSpanX operator()(DBlockSpanX ex, DBlockViewXVx fdistribu) const override;
};
