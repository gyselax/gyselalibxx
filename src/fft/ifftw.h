#pragma once

#include <array>
#include <complex>
#include <type_traits>

#include <blockview.h>
#include <fftw3.h>
#include <geometry.h>
#include <mdomain.h>
#include <non_uniform_mesh.h>

#include "ifft.h"

template <class... Tags>
class FftwInverseFourierTransform : public IInverseFourierTransform<Tags...>
{
public:
    FftwInverseFourierTransform() = default;

    FftwInverseFourierTransform(const FftwInverseFourierTransform& x) = default;

    FftwInverseFourierTransform(FftwInverseFourierTransform&& x) = default;

    ~FftwInverseFourierTransform() override = default;

    FftwInverseFourierTransform& operator=(const FftwInverseFourierTransform& x) = default;

    FftwInverseFourierTransform& operator=(FftwInverseFourierTransform&& x) = default;

    void operator()(
            BlockView<ProductMDomain<UniformMesh<Tags>...>, std::complex<double>> const& out_values,
            BlockView<ProductMDomain<NonUniformMesh<Fourier<Tags>>...>, std::complex<double>> const&
                    in_values) const noexcept override
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
                FFTW_BACKWARD,
                FFTW_ESTIMATE);
        fftw_execute(plan);

        fftw_destroy_plan(plan);
    }
};
