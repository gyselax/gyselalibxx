#pragma once

#include <complex>
#include <vector>

#include <blockview.h>
#include <geometry.h>
#include <mdomain.h>

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

    virtual std::vector<double> ifftshift(
            UniformMDomain<Fourier<Tags>...> const& dom_fx) const noexcept = 0;

    virtual void operator()(
            BlockView<UniformMDomain<Fourier<Tags>...>, std::complex<double>> const& out_values,
            BlockView<UniformMDomain<Tags...>, std::complex<double>> const& in_values)
            const noexcept = 0;
};
