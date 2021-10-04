#pragma once

#include <array>
#include <complex>
#include <type_traits>

#include <ddc/BlockSpan>
#include <ddc/NonUniformMesh>

#include <fftw3.h>
#include <geometry.h>

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
            BlockSpan<
                    std::complex<double>,
                    ProductMDomain<UniformMesh<Tags>...>,
                    std::experimental::layout_right> const& out_values,
            BlockSpan<
                    std::complex<double>,
                    ProductMDomain<NonUniformMesh<Fourier<Tags>>...>,
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
                FFTW_BACKWARD,
                FFTW_ESTIMATE);
        fftw_execute(plan);

        fftw_destroy_plan(plan);
    }
};
