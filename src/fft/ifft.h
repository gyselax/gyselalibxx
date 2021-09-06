#pragma once

#include <complex>

#include <blockview.h>
#include <geometry.h>
#include <mdomain.h>
#include <non_uniform_mesh.h>

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
            BlockView<
                    ProductMDomain<UniformMesh<Tags>...>,
                    std::complex<double>,
                    std::experimental::layout_right> const& out_values,
            BlockView<
                    ProductMDomain<NonUniformMesh<Fourier<Tags>>...>,
                    std::complex<double>,
                    std::experimental::layout_right> const& in_values) const noexcept = 0;
};
