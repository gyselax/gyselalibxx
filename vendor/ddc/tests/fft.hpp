// SPDX-License-Identifier: MIT

#include <experimental/mdspan>

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

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

// is_complex : trait to determine if type is Kokkos::complex<something>
template <typename T>
struct is_complex : std::false_type
{
};

template <typename T>
struct is_complex<Kokkos::complex<T>> : std::true_type
{
};

// LastSelector: returns a if Dim==Last, else b
template <typename T, typename Dim, typename Last>
constexpr T LastSelector(const T a, const T b)
{
    return std::is_same<Dim, Last>::value ? a : b;
}

template <typename T, typename Dim, typename First, typename Second, typename... Tail>
constexpr T LastSelector(const T a, const T b)
{
    return LastSelector<T, Dim, Second, Tail...>(a, b);
}

#if fftw_omp_AVAIL
template <typename MemorySpace, typename T>
using Allocator = typename std::conditional<
        std::is_same_v<
                MemorySpace,
                Kokkos::Serial::
                        memory_space> || std::is_same_v<MemorySpace, Kokkos::OpenMP::memory_space>,
        ddc::HostAllocator<T>,
        ddc::DeviceAllocator<T>>::type;

template <typename ExecSpace>
constexpr auto policy = [] {
    if constexpr (std::is_same<ExecSpace, Kokkos::Serial>::value) {
        return ddc::policies::serial_host;
    } else if constexpr (std::is_same<ExecSpace, Kokkos::OpenMP>::value) {
        return ddc::policies::parallel_host;
    } else {
        return ddc::policies::parallel_device;
    }
};
#else
template <typename MemorySpace, typename T>
using Allocator = typename std::conditional<
        std::is_same_v<MemorySpace, Kokkos::Serial::memory_space>,
        ddc::HostAllocator<T>,
        ddc::DeviceAllocator<T>>::type;

template <typename ExecSpace>
constexpr auto policy = [] {
    if constexpr (std::is_same<ExecSpace, Kokkos::Serial>::value) {
        return ddc::policies::serial_host;
    } else {
        return ddc::policies::parallel_device;
    }
};
#endif

// TODO:
// - FFT multidim but according to a subset of dimensions
template <typename ExecSpace, typename MemorySpace, typename Tin, typename Tout, typename... X>
static void test_fft()
{
    const double a = -10;
    const double b = 10;
    const std::size_t Nx = 64; // Optimal value is (b-a)^2/(2*pi)

    DDom<DDim<X>...> const x_mesh = DDom<DDim<X>...>(
            ddc::init_discrete_space(DDim<X>::
                                             init(ddc::Coordinate<X>(a + (b - a) / Nx / 2),
                                                  ddc::Coordinate<X>(b - (b - a) / Nx / 2),
                                                  DVect<DDim<X>>(Nx)))...);
    ddc::Chunk _f = ddc::Chunk(x_mesh, Allocator<MemorySpace, Tin>());
    ddc::ChunkSpan f = _f.span_view();
    ddc::for_each(
            policy<ExecSpace>(),
            ddc::get_domain<DDim<X>...>(f),
            DDC_LAMBDA(DElem<DDim<X>...> const e) {
                // f(e) = (Kokkos::cos(4*coordinate(ddc::select<DDim<X>>(e)))*...);
                // f(e) = ((Kokkos::sin(coordinate(ddc::select<DDim<X>>(e))+1e-20)/(coordinate(ddc::select<DDim<X>>(e))+1e-20))*...);
                f(e) = static_cast<Tin>(Kokkos::exp(
                        -(Kokkos::pow(coordinate(ddc::select<DDim<X>>(e)), 2) + ...) / 2));
            });

    ddc::init_fourier_space<X...>(x_mesh);
    DDom<DFDim<ddc::Fourier<X>>...> const k_mesh
            = ddc::FourierMesh(x_mesh, is_complex<Tin>::value && is_complex<Tout>::value);

    ddc::Chunk _f_bis = ddc::Chunk(ddc::get_domain<DDim<X>...>(f), Allocator<MemorySpace, Tin>());
    ddc::ChunkSpan f_bis = _f_bis.span_view();
    ddc::deepcopy(f_bis, f);

    ddc::Chunk _Ff = ddc::Chunk(k_mesh, Allocator<MemorySpace, Tout>());
    ddc::ChunkSpan Ff = _Ff.span_view();
    ddc::fft(ExecSpace(), Ff, f_bis, {ddc::FFT_Normalization::FULL});
    Kokkos::fence();

    // deepcopy of Ff because FFT C2R overwrites the input
    ddc::Chunk _Ff_bis = ddc::
            Chunk(ddc::get_domain<DFDim<ddc::Fourier<X>>...>(Ff), Allocator<MemorySpace, Tout>());
    ddc::ChunkSpan Ff_bis = _Ff_bis.span_view();
    ddc::deepcopy(Ff_bis, Ff);

    ddc::Chunk _FFf = ddc::Chunk(x_mesh, Allocator<MemorySpace, Tin>());
    ddc::ChunkSpan FFf = _FFf.span_view();
    ddc::ifft(ExecSpace(), FFf, Ff_bis, {ddc::FFT_Normalization::FULL});

    ddc::Chunk _f_host = ddc::Chunk(ddc::get_domain<DDim<X>...>(f), ddc::HostAllocator<Tin>());
    ddc::ChunkSpan f_host = _f_host.span_view();
    ddc::deepcopy(f_host, f);
#if 0
    std::cout << "\n input:\n";
    ddc::for_each(
            ddc::policies::serial_host,
            ddc::get_domain<DDim<X>...>(f_host),
            [=](DElem<DDim<X>...> const e) {
                (std::cout << ... << coordinate(ddc::select<DDim<X>>(e)))
                        << "->" << f_host(e) << ", ";
            });
#endif

    ddc::Chunk _Ff_host = ddc::
            Chunk(ddc::get_domain<DFDim<ddc::Fourier<X>>...>(Ff), ddc::HostAllocator<Tout>());
    ddc::ChunkSpan Ff_host = _Ff_host.span_view();
    ddc::deepcopy(Ff_host, Ff);
#if 0
    std::cout << "\n output:\n";
    ddc::for_each(
            ddc::policies::serial_host,
            ddc::get_domain<DFDim<ddc::Fourier<X>>...>(Ff_host),
            [=](DElem<DFDim<ddc::Fourier<X>>...> const e) {
                (std::cout << ... << coordinate(ddc::select<DFDim<ddc::Fourier<X>>>(e)))
                        << "->" << Kokkos::abs(Ff_host(e)) << " "
                        << Kokkos::exp(
                                   -(Kokkos::pow(coordinate(ddc::select<DFDim<ddc::Fourier<X>>>(e)), 2)
                                     + ...)
                                   / 2)
                        << ", ";
            });
#endif

    ddc::Chunk _FFf_host = ddc::Chunk(ddc::get_domain<DDim<X>...>(FFf), ddc::HostAllocator<Tin>());
    ddc::ChunkSpan FFf_host = _FFf_host.span_view();
    ddc::deepcopy(FFf_host, FFf);
#if 0
    std::cout << "\n iFFT(FFT):\n";
    ddc::for_each(
            ddc::policies::serial_host,
            ddc::get_domain<DDim<X>...>(FFf_host),
            [=](DElem<DDim<X>...> const e) {
                (std::cout << ... << coordinate(ddc::select<DDim<X>>(e)))
                        << "->" << Kokkos::abs(FFf_host(e)) << " " << Kokkos::abs(f_host(e))
                        << ", ";
            });
#endif

    double criterion = Kokkos::sqrt(ddc::transform_reduce(
            ddc::get_domain<DFDim<ddc::Fourier<X>>...>(Ff_host),
            0.,
            ddc::reducer::sum<double>(),
            [=](DElem<DFDim<ddc::Fourier<X>>...> const e) {
                return Kokkos::
                               pow(Kokkos::abs(Ff_host(e))
                                           - Kokkos::exp(
                                                   -(Kokkos::
                                                             pow(coordinate(ddc::select<
                                                                            DFDim<ddc::Fourier<X>>>(
                                                                         e)),
                                                                 2)
                                                     + ...)
                                                   / 2),
                                   2)
                       / (LastSelector<std::size_t, X, X...>(Nx / 2, Nx) * ...);
            }));

    double criterion2 = Kokkos::sqrt(ddc::transform_reduce(
            ddc::get_domain<DDim<X>...>(FFf_host),
            0.,
            ddc::reducer::sum<double>(),
            [=](DElem<DDim<X>...> const e) {
                return Kokkos::pow(Kokkos::abs(FFf_host(e)) - Kokkos::abs(f_host(e)), 2)
                       / Kokkos::pow(Nx, sizeof...(X));
            }));

    std::cout << "\n Distance between analytical prediction and numerical result : " << criterion;
    std::cout << "\n Distance between input and iFFT(FFT(input)) : " << criterion2;
    double epsilon = std::is_same_v<typename ddc::detail::fft::real_type<Tin>::type, double> ? 1e-15
                                                                                             : 1e-7;
    ASSERT_LE(criterion, epsilon);
    ASSERT_LE(criterion2, epsilon);
}
