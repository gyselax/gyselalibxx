// SPDX-License-Identifier: MIT

#pragma once

#include <complex>

#include <ddc/ddc.hpp>

#include <geometry.hpp>

template <class Tag>
class IInverseFourierTransform
{
public:
    virtual ~IInverseFourierTransform() = default;

    // Perform the normalized invFFT where the input is a complex and the output is a real
    virtual ChunkSpan<
            double,
            DiscreteDomain<UniformPointSampling<Tag>>,
            std::experimental::layout_right>
    operator()(
            ChunkSpan<
                    double,
                    DiscreteDomain<UniformPointSampling<Tag>>,
                    std::experimental::layout_right> out_values,
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformPointSampling<Fourier<Tag>>>,
                    std::experimental::layout_right> in_values) const noexcept = 0;

    // Perform the normalized invFFT where the input is a complex and the output is a complex
    virtual ChunkSpan<
            std::complex<double>,
            DiscreteDomain<UniformPointSampling<Tag>>,
            std::experimental::layout_right>
    operator()(
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<UniformPointSampling<Tag>>,
                    std::experimental::layout_right> out_values,
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformPointSampling<Fourier<Tag>>>,
                    std::experimental::layout_right> in_values) const noexcept = 0;
};
