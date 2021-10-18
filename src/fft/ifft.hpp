#pragma once

#include <complex>

#include <ddc/ChunkSpan>
#include <ddc/NonUniformDiscretization>

#include <geometry.hpp>

template <class Tag>
class IInverseFourierTransform
{
public:
    IInverseFourierTransform() = default;

    IInverseFourierTransform(const IInverseFourierTransform& x) = default;

    IInverseFourierTransform(IInverseFourierTransform&& x) = default;

    virtual ~IInverseFourierTransform() = default;

    IInverseFourierTransform& operator=(const IInverseFourierTransform& x) = default;

    IInverseFourierTransform& operator=(IInverseFourierTransform&& x) = default;

    // Perform the normalized invFFT where the input is a complex and the output is a real
    virtual void operator()(
            ChunkSpan<
                    double,
                    DiscreteDomain<UniformDiscretization<Tag>>,
                    std::experimental::layout_right> const& out_values,
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;

    // Perform the normalized invFFT where the input is a complex and the output is a complex
    virtual void operator()(
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<UniformDiscretization<Tag>>,
                    std::experimental::layout_right> const& out_values,
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;
};
