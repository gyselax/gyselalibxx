#pragma once

#include <complex>

#include <ddc/block_span.h>
#include <ddc/mdomain.h>
#include <ddc/non_uniform_mesh.h>

#include <geometry.h>

template <class... Tags>
class IInverseFourierTransform
{
public:
    IInverseFourierTransform() = default;

    IInverseFourierTransform(const IInverseFourierTransform& x) = default;

    IInverseFourierTransform(IInverseFourierTransform&& x) = default;

    virtual ~IInverseFourierTransform() = default;

    IInverseFourierTransform& operator=(const IInverseFourierTransform& x) = default;

    IInverseFourierTransform& operator=(IInverseFourierTransform&& x) = default;

    virtual void operator()(
            BlockSpan<
                    ProductMDomain<UniformMesh<Tags>...>,
                    std::complex<double>,
                    std::experimental::layout_right> const& out_values,
            BlockSpan<
                    ProductMDomain<NonUniformMesh<Fourier<Tags>>...>,
                    std::complex<double>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;
};
