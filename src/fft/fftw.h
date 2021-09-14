#pragma once

#include <array>
#include <complex>
#include <type_traits>

#include <ddc/BlockSpan>
#include <ddc/MDomain>
#include <ddc/NonUniformMesh>

#include <fftw3.h>
#include <geometry.h>

#include "fft.h"

template <class... Tags>
class FftwFourierTransform : public IFourierTransform<Tags...>
{
public:
    FftwFourierTransform() = default;

    FftwFourierTransform(const FftwFourierTransform& x) = default;

    FftwFourierTransform(FftwFourierTransform&& x) = default;

    ~FftwFourierTransform() override = default;

    FftwFourierTransform& operator=(const FftwFourierTransform& x) = default;

    FftwFourierTransform& operator=(FftwFourierTransform&& x) = default;

    NonUniformMesh<Fourier<Tags...>> compute_fourier_domain(
            ProductMDomain<UniformMesh<Tags>...> const& dom_x) const noexcept override
    {
        std::vector<double> freqs(dom_x.size());
        double const inv_Nd = 1. / (dom_x.size() * get<UniformMesh<Tags>...>(dom_x.mesh()).step());
        for (std::size_t ii = 0; ii <= dom_x.size() / 2; ++ii) {
            freqs[ii] = ii * inv_Nd;
        }
        for (std::size_t ii = dom_x.size() / 2 + 1; ii < dom_x.size(); ++ii) {
            freqs[ii] = -((dom_x.size() - ii) * inv_Nd);
        }
        return NonUniformMesh<Fourier<Tags>...>(freqs);
    }

    void operator()(
            BlockSpan<
                    ProductMDomain<NonUniformMesh<Fourier<Tags>>...>,
                    std::complex<double>,
                    std::experimental::layout_right> const& out_values,
            BlockSpan<
                    ProductMDomain<UniformMesh<Tags>...>,
                    std::complex<double>,
                    std::experimental::layout_right> const& in_values) const noexcept override
    {
        assert(in_values.extents().array() == out_values.extents().array());

        // It needs to be of type 'int'
        auto extents = out_values.extents();
        std::array<int, sizeof...(Tags)> n;
        for (std::size_t i = 0; i < extents.size(); ++i) {
            n[i] = extents[i];
        }

        fftw_plan plan = fftw_plan_dft(
                n.size(),
                n.data(),
                reinterpret_cast<fftw_complex*>(in_values.data()),
                reinterpret_cast<fftw_complex*>(out_values.data()),
                FFTW_FORWARD,
                FFTW_ESTIMATE);
        fftw_execute(plan);

        fftw_destroy_plan(plan);
    }
};
