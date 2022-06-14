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

    virtual typename NonUniformDiscretization<Fourier<Tag>>::template Impl<Kokkos::HostSpace>
    compute_fourier_domain(
            DiscreteDomain<UniformDiscretization<Tag>> const& dom_x) const noexcept = 0;

    // Perform FFT where the input is a real and the output is a complex
    virtual ChunkSpan<
            std::complex<double>,
            DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
            std::experimental::layout_right>
    operator()(
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
                    std::experimental::layout_right> out_values,
            ChunkSpan<
                    double,
                    DiscreteDomain<UniformDiscretization<Tag>>,
                    std::experimental::layout_right> in_values) const noexcept = 0;

    // Perform FFT where the input is a complex and output is a compleax
    virtual ChunkSpan<
            std::complex<double>,
            DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
            std::experimental::layout_right>
    operator()(
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
                    std::experimental::layout_right> out_values,
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<UniformDiscretization<Tag>>,
                    std::experimental::layout_right> in_values) const noexcept = 0;
};
