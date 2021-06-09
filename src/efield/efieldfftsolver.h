#pragma once

#include "fftw.h"
#include "iefieldsolver.h"
#include "ifftw.h"

class EfieldFftSolver : public IEfieldSolver
{
    IFourierTransform<Dim::X> const& m_fft = FftwFourierTransform<Dim::X>();

    IInverseFourierTransform<Dim::X> const& m_ifft = FftwInverseFourierTransform<Dim::X>();

public:
    DBlockSpanX operator()(DBlockSpanX ex, DBlockViewXVx fdistribu) const override;
};
