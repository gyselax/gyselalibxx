// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

#include <gtest/gtest.h>

template <typename X>
using DDim = ddc::UniformPointSampling<X>;

template <typename... DDim>
using DElem = ddc::DiscreteElement<DDim...>;

template <typename... DDim>
using DVect = ddc::DiscreteVector<DDim...>;

template <typename... DDim>
using DDom = ddc::DiscreteDomain<DDim...>;

template <typename Kx>
using DFDim = ddc::PeriodicSampling<Kx>;

template <typename ExecSpace>
constexpr auto policy = [] {
    if constexpr (std::is_same_v<ExecSpace, Kokkos::Serial>) {
        return ddc::policies::serial_host;
    }
#if fftw_omp_AVAIL
    else if constexpr (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
        return ddc::policies::parallel_host;
    }
#endif
    else {
        return ddc::policies::parallel_device;
    }
};

// TODO:
// - FFT multidim but according to a subset of dimensions
template <typename ExecSpace, typename MemorySpace, typename Tin, typename Tout, typename... X>
static void test_fft()
{
    constexpr bool full_fft
            = ddc::detail::fft::is_complex_v<Tin> && ddc::detail::fft::is_complex_v<Tout>;
    const double a = -10;
    const double b = 10;
    const std::size_t Nx = 64; // Optimal value is (b-a)^2/(2*pi)

    DDom<DDim<X>...> const x_mesh(
            ddc::init_discrete_space(DDim<X>::
                                             init(ddc::Coordinate<X>(a + (b - a) / Nx / 2),
                                                  ddc::Coordinate<X>(b - (b - a) / Nx / 2),
                                                  DVect<DDim<X>>(Nx)))...);
    ddc::init_fourier_space<X...>(x_mesh);
    DDom<DFDim<ddc::Fourier<X>>...> const k_mesh = ddc::FourierMesh(x_mesh, full_fft);

    ddc::Chunk _f(x_mesh, ddc::KokkosAllocator<Tin, MemorySpace>());
    ddc::ChunkSpan f = _f.span_view();
    ddc::for_each(
            policy<ExecSpace>(),
            f.domain(),
            DDC_LAMBDA(DElem<DDim<X>...> const e) {
                ddc::Real const xn2
                        = (Kokkos::pow(ddc::coordinate(ddc::select<DDim<X>>(e)), 2) + ...);
                f(e) = Kokkos::exp(-xn2 / 2);
            });
    ddc::Chunk f_bis_alloc(f.domain(), ddc::KokkosAllocator<Tin, MemorySpace>());
    ddc::ChunkSpan f_bis = f_bis_alloc.span_view();
    ddc::deepcopy(f_bis, f);

    ddc::Chunk Ff_alloc(k_mesh, ddc::KokkosAllocator<Tout, MemorySpace>());
    ddc::ChunkSpan Ff = Ff_alloc.span_view();
    ddc::fft(ExecSpace(), Ff, f_bis, {ddc::FFT_Normalization::FULL});
    Kokkos::fence();

    // deepcopy of Ff because FFT C2R overwrites the input
    ddc::Chunk Ff_bis_alloc(Ff.domain(), ddc::KokkosAllocator<Tout, MemorySpace>());
    ddc::ChunkSpan Ff_bis = Ff_bis_alloc.span_view();
    ddc::deepcopy(Ff_bis, Ff);

    ddc::Chunk FFf_alloc(f.domain(), ddc::KokkosAllocator<Tin, MemorySpace>());
    ddc::ChunkSpan FFf = FFf_alloc.span_view();
    ddc::ifft(ExecSpace(), FFf, Ff_bis, {ddc::FFT_Normalization::FULL});

    ddc::Chunk f_host_alloc(f.domain(), ddc::HostAllocator<Tin>());
    ddc::ChunkSpan f_host = f_host_alloc.span_view();
    ddc::deepcopy(f_host, f);

    ddc::Chunk Ff_host_alloc(Ff.domain(), ddc::HostAllocator<Tout>());
    ddc::ChunkSpan Ff_host = Ff_host_alloc.span_view();
    ddc::deepcopy(Ff_host, Ff);

    ddc::Chunk FFf_host_alloc(FFf.domain(), ddc::HostAllocator<Tin>());
    ddc::ChunkSpan FFf_host = FFf_host_alloc.span_view();
    ddc::deepcopy(FFf_host, FFf);

    auto const pow2 = DDC_LAMBDA(double x)
    {
        return x * x;
    };

    double criterion = Kokkos::sqrt(ddc::transform_reduce(
            Ff_host.domain(),
            0.,
            ddc::reducer::sum<double>(),
            [=](DElem<DFDim<ddc::Fourier<X>>...> const e) {
                double const xn2
                        = (pow2(ddc::coordinate(ddc::select<DFDim<ddc::Fourier<X>>>(e))) + ...);
                double const diff = Kokkos::abs(Ff_host(e)) - Kokkos::exp(-xn2 / 2);
                std::size_t const denom
                        = (ddc::detail::fft::LastSelector<std::size_t, X, X...>(Nx / 2, Nx) * ...);
                return pow2(diff) / denom;
            }));

    double criterion2 = Kokkos::sqrt(ddc::transform_reduce(
            FFf_host.domain(),
            0.,
            ddc::reducer::sum<double>(),
            [=](DElem<DDim<X>...> const e) {
                double const diff = Kokkos::abs(FFf_host(e)) - Kokkos::abs(f_host(e));
                return pow2(diff) / Kokkos::pow(Nx, sizeof...(X));
            }));
    double epsilon = std::is_same_v<ddc::detail::fft::real_type_t<Tin>, double> ? 1e-15 : 1e-7;
    EXPECT_LE(criterion, epsilon)
            << "Distance between analytical prediction and numerical result : " << criterion;
    EXPECT_LE(criterion2, epsilon)
            << "Distance between input and iFFT(FFT(input)) : " << criterion2;
}

template <typename ExecSpace, typename MemorySpace, typename Tin, typename Tout, typename X>
static void test_fft_norm(ddc::FFT_Normalization const norm)
{
    constexpr bool full_fft
            = ddc::detail::fft::is_complex_v<Tin> && ddc::detail::fft::is_complex_v<Tout>;

    DDom<DDim<X>> const x_mesh = ddc::init_discrete_space(DDim<X>::
                                                                  init(ddc::Coordinate<X>(-1. / 4),
                                                                       ddc::Coordinate<X>(1. / 4),
                                                                       DVect<DDim<X>>(2)));
    ddc::init_fourier_space<X>(x_mesh);
    DDom<DFDim<ddc::Fourier<X>>> const k_mesh = ddc::FourierMesh(x_mesh, full_fft);

    ddc::Chunk f_alloc = ddc::Chunk(x_mesh, ddc::KokkosAllocator<Tin, MemorySpace>());
    ddc::ChunkSpan f = f_alloc.span_view();
    ddc::for_each(
            policy<ExecSpace>(),
            f.domain(),
            DDC_LAMBDA(DElem<DDim<X>> const e) { f(e) = static_cast<Tin>(1); });


    ddc::Chunk f_bis_alloc = ddc::Chunk(f.domain(), ddc::KokkosAllocator<Tin, MemorySpace>());
    ddc::ChunkSpan f_bis = f_bis_alloc.span_view();
    ddc::deepcopy(f_bis, f);

    ddc::Chunk Ff_alloc = ddc::Chunk(k_mesh, ddc::KokkosAllocator<Tout, MemorySpace>());
    ddc::ChunkSpan Ff = Ff_alloc.span_view();
    ddc::fft(ExecSpace(), Ff, f_bis, {norm});
    Kokkos::fence();

    // deepcopy of Ff because FFT C2R overwrites the input
    ddc::Chunk Ff_bis_alloc = ddc::Chunk(Ff.domain(), ddc::KokkosAllocator<Tout, MemorySpace>());
    ddc::ChunkSpan Ff_bis = Ff_bis_alloc.span_view();
    ddc::deepcopy(Ff_bis, Ff);

    ddc::Chunk FFf_alloc = ddc::Chunk(x_mesh, ddc::KokkosAllocator<Tin, MemorySpace>());
    ddc::ChunkSpan FFf = FFf_alloc.span_view();
    ddc::ifft(ExecSpace(), FFf, Ff_bis, {norm});

    double const f_sum = ddc::transform_reduce(f.domain(), 0., ddc::reducer::sum<double>(), f);

    double Ff0_expected;
    double FFf_expected;
    switch (norm) {
    case ddc::FFT_Normalization::OFF:
        Ff0_expected = f_sum;
        FFf_expected = f_sum;
        break;
    case ddc::FFT_Normalization::FORWARD:
        Ff0_expected = 1;
        FFf_expected = 1;
        break;
    case ddc::FFT_Normalization::BACKWARD:
        Ff0_expected = f_sum;
        FFf_expected = 1;
        break;
    case ddc::FFT_Normalization::ORTHO:
        Ff0_expected = Kokkos::sqrt(f_sum);
        FFf_expected = 1;
        break;
    case ddc::FFT_Normalization::FULL:
        Ff0_expected = 1 / Kokkos::sqrt(2 * Kokkos::numbers::pi);
        FFf_expected = 1;
        break;
    }

    double epsilon = 1e-6;
    EXPECT_NEAR(Kokkos::abs(Ff(Ff.domain().front())), Ff0_expected, epsilon);
    EXPECT_NEAR(FFf(FFf.domain().front()), FFf_expected, epsilon);
    EXPECT_NEAR(FFf(FFf.domain().back()), FFf_expected, epsilon);
}
