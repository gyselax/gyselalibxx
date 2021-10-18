#pragma once

#include <array>
#include <complex>
#include <type_traits>

#include <ddc/ChunkSpan>
#include <ddc/NonUniformDiscretization>

#include <fftw3.h>
#include <geometry.hpp>

#include "ifft.hpp"

template <class Tag>
class FftwInverseFourierTransform : public IInverseFourierTransform<Tag>
{
public:
    FftwInverseFourierTransform() = default;

    FftwInverseFourierTransform(const FftwInverseFourierTransform& x) = default;

    FftwInverseFourierTransform(FftwInverseFourierTransform&& x) = default;

    ~FftwInverseFourierTransform() override = default;

    FftwInverseFourierTransform& operator=(const FftwInverseFourierTransform& x) = default;

    FftwInverseFourierTransform& operator=(FftwInverseFourierTransform&& x) = default;

    // Perform the normalized invFFT where the input is a complex and the output is a real
    void operator()(
            ChunkSpan<
                    double,
                    DiscreteDomain<UniformDiscretization<Tag>>,
                    std::experimental::layout_right> const& out_values,
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
                    std::experimental::layout_right> const& in_values) const noexcept override
    {
        assert(in_values.extents().array() == out_values.extents().array());

        // It needs to be of type 'int'
        auto extents = out_values.extents();
        std::array<int, 1> n;
        for (std::size_t i = 0; i < extents.size(); ++i) {
            n[i] = extents[i];
        }

        fftw_plan plan = fftw_plan_dft_c2r(
                n.size(),
                n.data(),
                reinterpret_cast<fftw_complex*>(in_values.data()),
                reinterpret_cast<double*>(out_values.data()),
                FFTW_ESTIMATE);
        fftw_execute(plan);

        fftw_destroy_plan(plan);

        // Normalization
        auto const out_domain = out_values.domain();
        for (auto i : out_domain) {
            out_values(i) = out_values(i) / out_domain.size();
        }
    }

    // Perform the normalized invFFT where the input is a complex and the output is a complex
    void operator()(
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<UniformDiscretization<Tag>>,
                    std::experimental::layout_right> const& out_values,
            ChunkSpan<
                    std::complex<double>,
                    DiscreteDomain<NonUniformDiscretization<Fourier<Tag>>>,
                    std::experimental::layout_right> const& in_values) const noexcept override
    {
        assert(in_values.extents().array() == out_values.extents().array());

        // It needs to be of type 'int'
        auto extents = out_values.extents();
        std::array<int, 1> n;
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

        // Normalization
        auto const out_domain = out_values.domain();
        for (auto i : out_domain) {
            out_values(i) /= out_domain.size();
        }
    }
};
