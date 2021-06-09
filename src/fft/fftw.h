#pragma once

#include <array>
#include <complex>
#include <type_traits>

#include <blockview.h>
#include <fftw3.h>
#include <geometry.h>
#include <mdomain.h>

#include "fft.h"

template <class... Tags>
class FftwFourierTransform : public IFourierTransform<Tags...>
{
public:
    FftwFourierTransform() = default;

    FftwFourierTransform(const FftwFourierTransform& x) = default;

    FftwFourierTransform(FftwFourierTransform&& x) = default;

    ~FftwFourierTransform() override = default;

    FftwFourierTransform& operator=(const FftwFourierTransform& x) = default;

    FftwFourierTransform& operator=(FftwFourierTransform&& x) = default;

    std::vector<double> ifftshift(
            UniformMDomain<Fourier<Tags>...> const& dom_fx) const noexcept override
    {
        std::vector<double> freqs(dom_fx.size());
        for (std::size_t ii = 0; ii <= dom_fx.size() / 2; ++ii) {
            freqs[ii] = dom_fx.to_real(ii);
        }
        for (std::size_t ii = dom_fx.size() / 2 + 1; ii < dom_fx.size(); ++ii) {
            freqs[ii] = dom_fx.to_real(ii) - dom_fx.rmax();
        }
        return freqs;
    }

    void operator()(
            BlockView<UniformMDomain<Fourier<Tags>...>, std::complex<double>> const& out_values,
            BlockView<UniformMDomain<Tags...>, std::complex<double>> const& in_values)
            const noexcept override
    {
        // shall be right layout
        std::array<int, out_values.rank()> n;
        for (std::size_t r = 0; r < out_values.rank(); ++r) {
            n[r] = out_values.extent(r);
        }

        fftw_plan plan = fftw_plan_dft(
                n.size(),
                n.data(),
                reinterpret_cast<fftw_complex*>(in_values.raw_view().data()),
                reinterpret_cast<fftw_complex*>(out_values.raw_view().data()),
                FFTW_FORWARD,
                FFTW_ESTIMATE);
        fftw_execute(plan);

        fftw_destroy_plan(plan);
    }
};
