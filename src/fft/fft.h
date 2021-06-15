#pragma once

#include <complex>
#include <vector>

#include <blockview.h>
#include <geometry.h>
#include <mdomain.h>
#include <nonuniformmesh.h>

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

    virtual NonUniformMDomain<Fourier<Tags>...> compute_fourier_domain(
            UniformMDomain<Tags...> const& dom_x) const noexcept = 0;

    virtual void operator()(
            BlockView<NonUniformMDomain<Fourier<Tags>...>, std::complex<double>> const& out_values,
            BlockView<UniformMDomain<Tags...>, std::complex<double>> const& in_values)
            const noexcept = 0;
};
