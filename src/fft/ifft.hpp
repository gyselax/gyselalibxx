// SPDX-License-Identifier: MIT

#pragma once

#include <complex>

#include <ddc/ddc.hpp>

#include "fft_tag.hpp"

template <class Tag>
class IInverseFourierTransform
{
public:
    virtual ~IInverseFourierTransform() = default;

    // Perform the normalized invFFT where the input is a complex and the output is a real
    virtual ddc::ChunkSpan<
            double,
            ddc::DiscreteDomain<ddc::UniformPointSampling<Tag>>,
            std::experimental::layout_right>
    operator()(
            ddc::ChunkSpan<
                    double,
                    ddc::DiscreteDomain<ddc::UniformPointSampling<Tag>>,
                    std::experimental::layout_right> out_values,
            ddc::ChunkSpan<
                    std::complex<double>,
                    ddc::DiscreteDomain<ddc::NonUniformPointSampling<Fourier<Tag>>>,
                    std::experimental::layout_right> in_values) const noexcept = 0;

    // Perform the normalized invFFT where the input is a complex and the output is a complex
    virtual ddc::ChunkSpan<
            std::complex<double>,
            ddc::DiscreteDomain<ddc::UniformPointSampling<Tag>>,
            std::experimental::layout_right>
    operator()(
            ddc::ChunkSpan<
                    std::complex<double>,
                    ddc::DiscreteDomain<ddc::UniformPointSampling<Tag>>,
                    std::experimental::layout_right> out_values,
            ddc::ChunkSpan<
                    std::complex<double>,
                    ddc::DiscreteDomain<ddc::NonUniformPointSampling<Fourier<Tag>>>,
                    std::experimental::layout_right> in_values) const noexcept = 0;
};
