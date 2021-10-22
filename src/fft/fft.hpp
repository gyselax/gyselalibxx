#pragma once

#include <complex>
#include <vector>

#include <ddc/ChunkSpan>
#include <ddc/NonUniformDiscretization>

#include <geometry.hpp>

template <class Tag>
class IFourierTransform
{
public:
    IFourierTransform() = default;

    IFourierTransform(IFourierTransform const& x) = default;

    IFourierTransform(IFourierTransform&& x) = default;

    virtual ~IFourierTransform() = default;

    IFourierTransform& operator=(IFourierTransform const& x) = default;

    IFourierTransform& operator=(IFourierTransform&& x) = default;

    virtual NonUniformDiscretization<Fourier<Tag>> compute_fourier_domain(
            DiscreteDomain<UniformDiscretization<Tag>> const& dom_x) const noexcept = 0;

    // Perform FFT where the input is a real and the output is a complex
    virtual void operator()(
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
                    std::experimental::layout_right> const& out_values,
            ChunkSpan<
                    double,
                    DiscreteDomain<UniformDiscretization<Tag>>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;

    // Perform FFT where the input is a complex and output is a compleax
    virtual void operator()(
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
                    std::experimental::layout_right> const& out_values,
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<UniformDiscretization<Tag>>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;
};
