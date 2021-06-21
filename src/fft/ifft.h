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
            BlockView<UniformMDomain<Tags...>, std::complex<double>> const& out_values,
            BlockView<NonUniformMDomain<Fourier<Tags>...>, std::complex<double>> const& in_values)
            const noexcept = 0;
};
