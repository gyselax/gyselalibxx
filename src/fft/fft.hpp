// SPDX-License-Identifier: MIT

#pragma once

#include <complex>
#include <vector>

#include <ddc/ddc.hpp>

#include <geometry.hpp>

template <class Tag>
class IFourierTransform
{
public:
    virtual ~IFourierTransform() = default;

    virtual typename ddc::NonUniformPointSampling<Fourier<Tag>>::template Impl<Kokkos::HostSpace>
    compute_fourier_domain(
            ddc::DiscreteDomain<ddc::UniformPointSampling<Tag>> const& dom_x) const noexcept = 0;

    // Perform FFT where the input is a real and the output is a complex
    virtual ddc::ChunkSpan<
            std::complex<double>,
            ddc::DiscreteDomain<ddc::NonUniformPointSampling<Fourier<Tag>>>,
            std::experimental::layout_right>
    operator()(
            ddc::ChunkSpan<
                    std::complex<double>,
                    ddc::DiscreteDomain<ddc::NonUniformPointSampling<Fourier<Tag>>>,
                    std::experimental::layout_right> out_values,
            ddc::ChunkSpan<
                    double,
                    ddc::DiscreteDomain<ddc::UniformPointSampling<Tag>>,
                    std::experimental::layout_right> in_values) const noexcept = 0;

    // Perform FFT where the input is a complex and output is a compleax
    virtual ddc::ChunkSpan<
            std::complex<double>,
            ddc::DiscreteDomain<ddc::NonUniformPointSampling<Fourier<Tag>>>,
            std::experimental::layout_right>
    operator()(
            ddc::ChunkSpan<
                    std::complex<double>,
                    ddc::DiscreteDomain<ddc::NonUniformPointSampling<Fourier<Tag>>>,
                    std::experimental::layout_right> out_values,
            ddc::ChunkSpan<
                    std::complex<double>,
                    ddc::DiscreteDomain<ddc::UniformPointSampling<Tag>>,
                    std::experimental::layout_right> in_values) const noexcept = 0;
};
