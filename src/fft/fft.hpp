#pragma once

#include <complex>
#include <vector>

#include <ddc/BlockSpan>
#include <ddc/NonUniformMesh>

#include <geometry.hpp>

template <class Tag>
class IFourierTransform
{
public:
    IFourierTransform() = default;

    IFourierTransform(const IFourierTransform& x) = default;

    IFourierTransform(IFourierTransform&& x) = default;

    virtual ~IFourierTransform() = default;

    IFourierTransform& operator=(const IFourierTransform& x) = default;

    IFourierTransform& operator=(IFourierTransform&& x) = default;

    virtual NonUniformMesh<Fourier<Tag>> compute_fourier_domain(
            ProductMDomain<UniformMesh<Tag>> const& dom_x) const noexcept = 0;

    // Perform FFT where the input is a real and the output is a complex
    virtual void operator()(
            BlockSpan<
                    std::complex<double>,
                    ProductMDomain<NonUniformMesh<Fourier<Tag>>>,
                    std::experimental::layout_right> const& out_values,
            BlockSpan<
                    double,
                    ProductMDomain<UniformMesh<Tag>>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;

    // Perform FFT where the input is a complex and output is a compleax
    virtual void operator()(
            BlockSpan<
                    std::complex<double>,
                    ProductMDomain<NonUniformMesh<Fourier<Tag>>>,
                    std::experimental::layout_right> const& out_values,
            BlockSpan<
                    std::complex<double>,
                    ProductMDomain<UniformMesh<Tag>>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;
};
