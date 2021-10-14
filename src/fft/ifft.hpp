#pragma once

#include <complex>

#include <ddc/BlockSpan>
#include <ddc/NonUniformMesh>

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
            BlockSpan<
                    double,
                    ProductMDomain<UniformMesh<Tag>>,
                    std::experimental::layout_right> const& out_values,
            BlockSpan<
                    std::complex<double>,
                    ProductMDomain<NonUniformMesh<Fourier<Tag>>>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;

    // Perform the normalized invFFT where the input is a complex and the output is a complex
    virtual void operator()(
            BlockSpan<
                    std::complex<double>,
                    ProductMDomain<UniformMesh<Tag>>,
                    std::experimental::layout_right> const& out_values,
            BlockSpan<
                    std::complex<double>,
                    ProductMDomain<NonUniformMesh<Fourier<Tag>>>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;
};
