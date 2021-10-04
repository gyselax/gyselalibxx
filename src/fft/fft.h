#pragma once

#include <complex>
#include <vector>

#include <ddc/BlockSpan>
#include <ddc/NonUniformMesh>

#include <geometry.h>

template <class... Tags>
class IFourierTransform
{
public:
    IFourierTransform() = default;

    IFourierTransform(const IFourierTransform& x) = default;

    IFourierTransform(IFourierTransform&& x) = default;

    virtual ~IFourierTransform() = default;

    IFourierTransform& operator=(const IFourierTransform& x) = default;

    IFourierTransform& operator=(IFourierTransform&& x) = default;

    virtual NonUniformMesh<Fourier<Tags...>> compute_fourier_domain(
            ProductMDomain<UniformMesh<Tags>...> const& dom_x) const noexcept = 0;

    virtual void operator()(
            BlockSpan<
                    std::complex<double>,
                    ProductMDomain<NonUniformMesh<Fourier<Tags>>...>,
                    std::experimental::layout_right> const& out_values,
            BlockSpan<
                    std::complex<double>,
                    ProductMDomain<UniformMesh<Tags>...>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;
};
