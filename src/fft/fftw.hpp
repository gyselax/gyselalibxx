// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <complex>
#include <type_traits>

#include <ddc/ddc.hpp>

#include <fftw3.h>
#include <geometry.hpp>

#include "fft.hpp"

template <class Tag>
class FftwFourierTransform : public IFourierTransform<Tag>
{
public:
    FftwFourierTransform() = default;

    ~FftwFourierTransform() override = default;

    NonUniformDiscretization<Fourier<Tag>> compute_fourier_domain(
            DiscreteDomain<UniformDiscretization<Tag>> const& dom_x) const noexcept override
    {
        std::vector<double> freqs(dom_x.size());
        double const inv_Nd = 1. / (dom_x.size() * step<UniformDiscretization<Tag>>());
        for (std::size_t ii = 0; ii <= dom_x.size() / 2; ++ii) {
            freqs[ii] = ii * inv_Nd;
        }
        for (std::size_t ii = dom_x.size() / 2 + 1; ii < dom_x.size(); ++ii) {
            freqs[ii] = -((dom_x.size() - ii) * inv_Nd);
        }
        return NonUniformDiscretization<Fourier<Tag>>(freqs);
    }

    // Perform FFT where the input is a real and the output is a complex
    ChunkSpan<
            std::complex<double>,
            DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
            std::experimental::layout_right>
    operator()(
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
                    std::experimental::layout_right> const out_values,
            ChunkSpan<
                    double,
                    DiscreteDomain<UniformDiscretization<Tag>>,
                    std::experimental::layout_right> const in_values) const noexcept override
    {
        assert(in_values.extents().array() == out_values.extents().array());

        // It needs to be of type 'int'
        auto extents = out_values.extents();
        std::array<int, 1> n;
        for (std::size_t i = 0; i < extents.size(); ++i) {
            n[i] = extents[i];
        }

        fftw_plan plan = fftw_plan_dft_r2c(
                n.size(),
                n.data(),
                reinterpret_cast<double*>(in_values.data()),
                reinterpret_cast<fftw_complex*>(out_values.data()),
                FFTW_ESTIMATE);
        fftw_execute(plan);

        fftw_destroy_plan(plan);

        return out_values;
    }

    // Perform FFT where the input is a complex and the output is a complex
    ChunkSpan<
            std::complex<double>,
            DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
            std::experimental::layout_right>
    operator()(
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
                    std::experimental::layout_right> const out_values,
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<UniformDiscretization<Tag>>,
                    std::experimental::layout_right> const in_values) const noexcept override
    {
        assert(in_values.extents().array() == out_values.extents().array());

        // It needs to be of type 'int'
        auto const extents = out_values.extents();
        std::array<int, 1> n;
        for (std::size_t i = 0; i < extents.size(); ++i) {
            n[i] = extents[i];
        }

        fftw_plan plan = fftw_plan_dft(
                n.size(),
                n.data(),
                reinterpret_cast<fftw_complex*>(in_values.data()),
                reinterpret_cast<fftw_complex*>(out_values.data()),
                FFTW_FORWARD,
                FFTW_ESTIMATE);
        fftw_execute(plan);

        fftw_destroy_plan(plan);

        return out_values;
    }
};
