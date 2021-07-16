#pragma once

#include <array>
#include <complex>
#include <type_traits>

#include <blockview.h>
#include <fftw3.h>
#include <geometry.h>
#include <mdomain.h>
#include <non_uniform_mesh.h>

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
            BlockView<ProductMDomain<NonUniformMesh<Fourier<Tags>>...>, std::complex<double>> const&
                    out_values,
            BlockView<ProductMDomain<UniformMesh<Tags>...>, std::complex<double>> const& in_values)
            const noexcept override
    {
        assert(in_values.is_contiguous());
        assert(in_values.is_unique());
        assert(out_values.is_contiguous());
        assert(out_values.is_unique());

        // shall be right layout
        std::array<int, out_values.rank()> n;
        for (std::size_t r = 0; r < out_values.rank(); ++r) {
            n[r] = out_values.extent(r);
        }

        fftw_plan plan = fftw_plan_dft(
                n.size(),
                n.data(),
                reinterpret_cast<fftw_complex*>(in_values.raw_view().data()),
                reinterpret_cast<fftw_complex*>(out_values.raw_view().data()),
                FFTW_FORWARD,
                FFTW_ESTIMATE);
        fftw_execute(plan);

        fftw_destroy_plan(plan);
    }
};
