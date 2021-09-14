#pragma once

#include <complex>
#include <vector>

#include <block_span.h>
#include <geometry.h>
#include <mdomain.h>
#include <non_uniform_mesh.h>

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
                    ProductMDomain<NonUniformMesh<Fourier<Tags>>...>,
                    std::complex<double>,
                    std::experimental::layout_right> const& out_values,
            BlockSpan<
                    ProductMDomain<UniformMesh<Tags>...>,
                    std::complex<double>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;
};
