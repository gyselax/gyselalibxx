// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <complex>
#include <type_traits>

#include <ddc/ddc.hpp>

#include <fftw3.h>
#include <geometry.hpp>

#include "ifft.hpp"

template <class Tag>
class FftwInverseFourierTransform : public IInverseFourierTransform<Tag>
{
public:
    FftwInverseFourierTransform() = default;

    ~FftwInverseFourierTransform() override = default;

    // Perform the normalized invFFT where the input is a complex and the output is a real
    ddc::ChunkSpan<
            double,
            ddc::DiscreteDomain<ddc::UniformPointSampling<Tag>>,
            std::experimental::layout_right>
    operator()(
            ddc::ChunkSpan<
                    double,
                    ddc::DiscreteDomain<ddc::UniformPointSampling<Tag>>,
                    std::experimental::layout_right> const out_values,
            ddc::ChunkSpan<
                    std::complex<double>,
                    ddc::DiscreteDomain<ddc::NonUniformPointSampling<Fourier<Tag>>>,
                    std::experimental::layout_right> const in_values) const noexcept override
    {
        assert(in_values.extents().value() == out_values.extents().value());

        // It needs to be of type 'int'
        int const n = out_values.extents().value();

        fftw_plan plan = fftw_plan_dft_c2r(
                1,
                &n,
                reinterpret_cast<fftw_complex*>(in_values.data()),
                reinterpret_cast<double*>(out_values.data()),
                FFTW_ESTIMATE);
        fftw_execute(plan);

        fftw_destroy_plan(plan);

        // Normalization
        auto const out_domain = out_values.domain();
        for (auto const i : out_domain) {
            out_values(i) = out_values(i) / out_domain.size();
        }

        return out_values;
    }

    // Perform the normalized invFFT where the input is a complex and the output is a complex
    ddc::ChunkSpan<
            std::complex<double>,
            ddc::DiscreteDomain<ddc::UniformPointSampling<Tag>>,
            std::experimental::layout_right>
    operator()(
            ddc::ChunkSpan<
                    std::complex<double>,
                    ddc::DiscreteDomain<ddc::UniformPointSampling<Tag>>,
                    std::experimental::layout_right> const out_values,
            ddc::ChunkSpan<
                    std::complex<double>,
                    ddc::DiscreteDomain<ddc::NonUniformPointSampling<Fourier<Tag>>>,
                    std::experimental::layout_right> const in_values) const noexcept override
    {
        assert(in_values.extents().value() == out_values.extents().value());

        // It needs to be of type 'int'
        int const n = out_values.extents().value();

        fftw_plan plan = fftw_plan_dft(
                1,
                &n,
                reinterpret_cast<fftw_complex*>(in_values.data()),
                reinterpret_cast<fftw_complex*>(out_values.data()),
                FFTW_BACKWARD,
                FFTW_ESTIMATE);
        fftw_execute(plan);

        fftw_destroy_plan(plan);

        // Normalization
        auto const out_domain = out_values.domain();
        for (auto const i : out_domain) {
            out_values(i) /= out_domain.size();
        }

        return out_values;
    }
};
