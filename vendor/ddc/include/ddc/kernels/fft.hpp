// SPDX-License-Identifier: MIT

#pragma once

#include <memory>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>

#if fftw_AVAIL
#include <fftw3.h>
#endif

#if cufft_AVAIL
#include <cufft.h>
#endif

#if hipfft_AVAIL
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hipfft/hipfft.h>
#endif

#if fftw_AVAIL
static_assert(sizeof(fftwf_complex) == sizeof(Kokkos::complex<float>));
static_assert(alignof(fftwf_complex) <= alignof(Kokkos::complex<float>));

static_assert(sizeof(fftw_complex) == sizeof(Kokkos::complex<double>));
static_assert(alignof(fftw_complex) <= alignof(Kokkos::complex<double>));

static_assert(sizeof(fftwl_complex) == sizeof(Kokkos::complex<long double>));
static_assert(alignof(fftwl_complex) <= alignof(Kokkos::complex<long double>));
#endif

#if cufft_AVAIL
static_assert(sizeof(cufftComplex) == sizeof(Kokkos::complex<float>));
static_assert(alignof(cufftComplex) <= alignof(Kokkos::complex<float>));

static_assert(sizeof(cufftDoubleComplex) == sizeof(Kokkos::complex<double>));
static_assert(alignof(cufftDoubleComplex) <= alignof(Kokkos::complex<double>));
#endif

#if hipfft_AVAIL
static_assert(sizeof(hipfftComplex) == sizeof(Kokkos::complex<float>));
static_assert(alignof(hipfftComplex) <= alignof(Kokkos::complex<float>));

static_assert(sizeof(hipfftDoubleComplex) == sizeof(Kokkos::complex<double>));
static_assert(alignof(hipfftDoubleComplex) <= alignof(Kokkos::complex<double>));
#endif

namespace ddc {
// TODO : maybe transfert this somewhere else because Fourier space is not specific to FFT
template <typename Dim>
struct Fourier;

// named arguments for FFT (and their default values)
enum class FFT_Direction { FORWARD, BACKWARD };
enum class FFT_Normalization { OFF, FORWARD, BACKWARD, ORTHO, FULL };
} // namespace ddc

namespace ddc::detail::fft {
template <typename T>
struct real_type
{
    using type = T;
};

template <typename T>
struct real_type<Kokkos::complex<T>>
{
    using type = T;
};

template <typename T>
using real_type_t = typename real_type<T>::type;

// is_complex : trait to determine if type is Kokkos::complex<something>
template <typename T>
struct is_complex : std::false_type
{
};

template <typename T>
struct is_complex<Kokkos::complex<T>> : std::true_type
{
};

template <typename T>
constexpr bool is_complex_v = is_complex<T>::value;

// LastSelector: returns a if Dim==Last, else b
template <typename T, typename Dim, typename Last>
KOKKOS_FUNCTION constexpr T LastSelector(const T a, const T b)
{
    return std::is_same_v<Dim, Last> ? a : b;
}

template <typename T, typename Dim, typename First, typename Second, typename... Tail>
KOKKOS_FUNCTION constexpr T LastSelector(const T a, const T b)
{
    return LastSelector<T, Dim, Second, Tail...>(a, b);
}

// transform_type : trait to determine the type of transformation (R2C, C2R, C2C...) <- no information about base type (float or double)
enum class TransformType { R2R, R2C, C2R, C2C };

template <typename T1, typename T2>
struct transform_type
{
    static constexpr TransformType value = TransformType::R2R;
};

template <typename T1, typename T2>
struct transform_type<T1, Kokkos::complex<T2>>
{
    static constexpr TransformType value = TransformType::R2C;
};

template <typename T1, typename T2>
struct transform_type<Kokkos::complex<T1>, T2>
{
    static constexpr TransformType value = TransformType::C2R;
};

template <typename T1, typename T2>
struct transform_type<Kokkos::complex<T1>, Kokkos::complex<T2>>
{
    static constexpr TransformType value = TransformType::C2C;
};

template <typename T1, typename T2>
constexpr TransformType transform_type_v = transform_type<T1, T2>::value;

#if fftw_AVAIL
// _fftw_type : compatible with both single and double precision
template <typename T>
struct _fftw_type
{
    using type = T;
};

template <typename T>
struct _fftw_type<Kokkos::complex<T>>
{
    using type = std::
            conditional_t<std::is_same_v<real_type_t<T>, float>, fftwf_complex, fftw_complex>;
};

// _fftw_plan : compatible with both single and double precision
template <typename T>
using _fftw_plan = std::conditional_t<std::is_same_v<real_type_t<T>, float>, fftwf_plan, fftw_plan>;

// _fftw_plan_many_dft : templated function working for all types of transformation
template <typename Tin, typename Tout, typename... Args, typename PenultArg, typename LastArg>
_fftw_plan<Tin> _fftw_plan_many_dft(
        [[maybe_unused]] PenultArg penultArg,
        LastArg lastArg,
        Args... args)
{ // Ugly, penultArg and lastArg are passed before the rest because of a limitation of C++ (parameter packs must be last arguments)
    const TransformType transformType = transform_type_v<Tin, Tout>;
    if constexpr (transformType == TransformType::R2C && std::is_same_v<Tin, float>)
        return fftwf_plan_many_dft_r2c(args..., lastArg);
    else if constexpr (transformType == TransformType::R2C && std::is_same_v<Tin, double>)
        return fftw_plan_many_dft_r2c(args..., lastArg);
    else if constexpr (transformType == TransformType::C2R && std::is_same_v<Tout, float>)
        return fftwf_plan_many_dft_c2r(args..., lastArg);
    else if constexpr (transformType == TransformType::C2R && std::is_same_v<Tout, double>)
        return fftw_plan_many_dft_c2r(args..., lastArg);
    else if constexpr (
            transformType == TransformType::C2C && std::is_same_v<Tin, Kokkos::complex<float>>)
        return fftwf_plan_many_dft(args..., penultArg, lastArg);
    else if constexpr (
            transformType == TransformType::C2C && std::is_same_v<Tin, Kokkos::complex<double>>)
        return fftw_plan_many_dft(args..., penultArg, lastArg);
    // else constexpr
    //   static_assert(false, "Transform type not supported");
}

#endif
#if cufft_AVAIL
// _cufft_type : compatible with both single and double precision
template <typename T>
struct _cufft_type
{
    using type = std::conditional_t<std::is_same_v<T, float>, cufftReal, cufftDoubleReal>;
};

template <typename T>
struct _cufft_type<Kokkos::complex<T>>
{
    using type = std::
            conditional_t<std::is_same_v<real_type_t<T>, float>, cufftComplex, cufftDoubleComplex>;
};

// cufft_transform_type : argument passed in the cufftMakePlan function
template <typename Tin, typename Tout>
constexpr auto cufft_transform_type()
{
    const TransformType transformType = transform_type_v<Tin, Tout>;
    if constexpr (transformType == TransformType::R2C && std::is_same_v<Tin, float>)
        return CUFFT_R2C;
    else if constexpr (transformType == TransformType::R2C && std::is_same_v<Tin, double>)
        return CUFFT_D2Z;
    else if constexpr (transformType == TransformType::C2R && std::is_same_v<Tout, float>)
        return CUFFT_C2R;
    else if constexpr (transformType == TransformType::C2R && std::is_same_v<Tout, double>)
        return CUFFT_Z2D;
    else if constexpr (
            transformType == TransformType::C2C && std::is_same_v<Tin, Kokkos::complex<float>>)
        return CUFFT_C2C;
    else if constexpr (
            transformType == TransformType::C2C && std::is_same_v<Tin, Kokkos::complex<double>>)
        return CUFFT_Z2Z;
    // else constexpr
    //	static_assert(false, "Transform type not supported");
}

// cufftExec : argument passed in the cufftMakePlan function
// _fftw_plan_many_dft : templated function working for all types of transformation
template <typename Tin, typename Tout, typename... Args, typename LastArg>
cufftResult _cufftExec([[maybe_unused]] LastArg lastArg, Args... args)
{ // Ugly for same reason as fftw
    const TransformType transformType = transform_type_v<Tin, Tout>;
    if constexpr (transformType == TransformType::R2C && std::is_same_v<Tin, float>)
        return cufftExecR2C(args...);
    else if constexpr (transformType == TransformType::R2C && std::is_same_v<Tin, double>)
        return cufftExecD2Z(args...);
    else if constexpr (transformType == TransformType::C2R && std::is_same_v<Tout, float>)
        return cufftExecC2R(args...);
    else if constexpr (transformType == TransformType::C2R && std::is_same_v<Tout, double>)
        return cufftExecZ2D(args...);
    else if constexpr (
            transformType == TransformType::C2C && std::is_same_v<Tin, Kokkos::complex<float>>)
        return cufftExecC2C(args..., lastArg);
    else if constexpr (
            transformType == TransformType::C2C && std::is_same_v<Tin, Kokkos::complex<double>>)
        return cufftExecZ2Z(args..., lastArg);
    // else constexpr
    //   static_assert(false, "Transform type not supported");
}
#endif
#if hipfft_AVAIL
// _hipfft_type : compatible with both single and double precision
template <typename T>
struct _hipfft_type
{
    using type = std::conditional_t<std::is_same_v<T, float>, hipfftReal, hipfftDoubleReal>;
};

template <typename T>
struct _hipfft_type<Kokkos::complex<T>>
{
    using type = std::conditional_t<
            std::is_same_v<real_type_t<T>, float>,
            hipfftComplex,
            hipfftDoubleComplex>;
};

// hipfft_transform_type : argument passed in the hipfftMakePlan function
template <typename Tin, typename Tout>
constexpr auto hipfft_transform_type()
{
    const TransformType transformType = transform_type_v<Tin, Tout>;
    if constexpr (transformType == TransformType::R2C && std::is_same_v<Tin, float>)
        return HIPFFT_R2C;
    else if constexpr (transformType == TransformType::R2C && std::is_same_v<Tin, double>)
        return HIPFFT_D2Z;
    else if constexpr (transformType == TransformType::C2R && std::is_same_v<Tout, float>)
        return HIPFFT_C2R;
    else if constexpr (transformType == TransformType::C2R && std::is_same_v<Tout, double>)
        return HIPFFT_Z2D;
    else if constexpr (
            transformType == TransformType::C2C && std::is_same_v<Tin, Kokkos::complex<float>>)
        return HIPFFT_C2C;
    else if constexpr (
            transformType == TransformType::C2C && std::is_same_v<Tin, Kokkos::complex<double>>)
        return HIPFFT_Z2Z;
    // else constexpr
    //	static_assert(false, "Transform type not supported");
}

// hipfftExec : argument passed in the hipfftMakePlan function
// _fftw_plan_many_dft : templated function working for all types of transformation
template <typename Tin, typename Tout, typename... Args, typename LastArg>
hipfftResult _hipfftExec([[maybe_unused]] LastArg lastArg, Args... args)
{
    const TransformType transformType = transform_type_v<Tin, Tout>;
    if constexpr (transformType == TransformType::R2C && std::is_same_v<Tin, float>)
        return hipfftExecR2C(args...);
    else if constexpr (transformType == TransformType::R2C && std::is_same_v<Tin, double>)
        return hipfftExecD2Z(args...);
    else if constexpr (transformType == TransformType::C2R && std::is_same_v<Tout, float>)
        return hipfftExecC2R(args...);
    else if constexpr (transformType == TransformType::C2R && std::is_same_v<Tout, double>)
        return hipfftExecZ2D(args...);
    else if constexpr (
            transformType == TransformType::C2C && std::is_same_v<Tin, Kokkos::complex<float>>)
        return hipfftExecC2C(args..., lastArg);
    else if constexpr (
            transformType == TransformType::C2C && std::is_same_v<Tin, Kokkos::complex<double>>)
        return hipfftExecZ2Z(args..., lastArg);
    // else constexpr
    //   static_assert(false, "Transform type not supported");
}
#endif

struct kwArgs_core
{
    ddc::FFT_Direction
            direction; // Only effective for C2C transform and for normalization BACKWARD and FORWARD
    ddc::FFT_Normalization normalization;
};

// N,a,b from x_mesh
template <typename Dim, typename... X>
int N(ddc::DiscreteDomain<ddc::UniformPointSampling<X>...> x_mesh)
{
    return ddc::get<ddc::UniformPointSampling<Dim>>(x_mesh.extents());
}

template <typename Dim, typename... X>
double a(ddc::DiscreteDomain<ddc::UniformPointSampling<X>...> x_mesh)
{
    return ((2 * N<Dim>(x_mesh) - 1)
                    * coordinate(ddc::select<ddc::UniformPointSampling<Dim>>(x_mesh).front())
            - coordinate(ddc::select<ddc::UniformPointSampling<Dim>>(x_mesh).back()))
           / 2 / (N<Dim>(x_mesh) - 1);
}

template <typename Dim, typename... X>
double b(ddc::DiscreteDomain<ddc::UniformPointSampling<X>...> x_mesh)
{
    return ((2 * N<Dim>(x_mesh) - 1)
                    * coordinate(ddc::select<ddc::UniformPointSampling<Dim>>(x_mesh).back())
            - coordinate(ddc::select<ddc::UniformPointSampling<Dim>>(x_mesh).front()))
           / 2 / (N<Dim>(x_mesh) - 1);
}

// core
template <typename Tin, typename Tout, typename ExecSpace, typename MemorySpace, typename... X>
void core(
        ExecSpace& execSpace,
        Tout* out_data,
        Tin* in_data,
        ddc::DiscreteDomain<ddc::UniformPointSampling<X>...> mesh,
        const kwArgs_core& kwargs)
{
    static_assert(
            std::is_same_v<real_type_t<Tin>, float> || std::is_same_v<real_type_t<Tin>, double>,
            "Base type of Tin (and Tout) must be float or double.");
    static_assert(
            std::is_same_v<real_type_t<Tin>, real_type_t<Tout>>,
            "Types Tin and Tout must be based on same type (float or double)");
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible,
            "MemorySpace has to be accessible for ExecutionSpace.");

    std::array<int, sizeof...(X)> n
            = {(int)ddc::get<ddc::UniformPointSampling<X>>(mesh.extents())...};
    int idist = 1;
    int odist = 1;
    for (size_t i = 0; i < sizeof...(X); i++) {
        idist = transform_type_v<Tin, Tout> == TransformType::C2R && i == sizeof...(X) - 1
                        ? idist * (n[i] / 2 + 1)
                        : idist * n[i];
        odist = transform_type_v<Tin, Tout> == TransformType::R2C && i == sizeof...(X) - 1
                        ? odist * (n[i] / 2 + 1)
                        : odist * n[i];
    }

    if constexpr (false) {
    } // Trick to get only else if
#if fftw_AVAIL
    else if constexpr (std::is_same_v<ExecSpace, Kokkos::Serial>) {
        _fftw_plan<Tin> plan = _fftw_plan_many_dft<Tin, Tout>(
                kwargs.direction == ddc::FFT_Direction::FORWARD ? FFTW_FORWARD : FFTW_BACKWARD,
                FFTW_ESTIMATE,
                (int)sizeof...(X),
                n.data(),
                1,
                reinterpret_cast<typename _fftw_type<Tin>::type*>(in_data),
                (int*)NULL,
                1,
                idist,
                reinterpret_cast<typename _fftw_type<Tout>::type*>(out_data),
                (int*)NULL,
                1,
                odist);
        if constexpr (std::is_same_v<real_type_t<Tin>, float>) {
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
        } else {
            fftw_execute(plan);
            fftw_destroy_plan(plan);
        }
        // std::cout << "performed with fftw";
    }
#endif
#if fftw_omp_AVAIL
    else if constexpr (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
        if constexpr (std::is_same_v<real_type_t<Tin>, float>) {
            fftwf_init_threads();
            fftwf_plan_with_nthreads(execSpace.concurrency());
        } else {
            fftw_init_threads();
            fftw_plan_with_nthreads(execSpace.concurrency());
        }
        _fftw_plan<Tin> plan = _fftw_plan_many_dft<Tin, Tout>(
                kwargs.direction == ddc::FFT_Direction::FORWARD ? FFTW_FORWARD : FFTW_BACKWARD,
                FFTW_ESTIMATE,
                (int)sizeof...(X),
                n.data(),
                1,
                reinterpret_cast<typename _fftw_type<Tin>::type*>(in_data),
                (int*)NULL,
                1,
                idist,
                reinterpret_cast<typename _fftw_type<Tout>::type*>(out_data),
                (int*)NULL,
                1,
                odist);
        if constexpr (std::is_same_v<real_type_t<Tin>, float>) {
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
        } else {
            fftw_execute(plan);
            fftw_destroy_plan(plan);
        }
        // std::cout << "performed with fftw_omp";
    }
#endif
#if cufft_AVAIL
    else if constexpr (std::is_same_v<ExecSpace, Kokkos::Cuda>) {
        cudaStream_t stream = execSpace.cuda_stream();

        cufftHandle unmanaged_plan = -1;
        cufftResult cufft_rt = cufftCreate(&unmanaged_plan);

        if (cufft_rt != CUFFT_SUCCESS)
            throw std::runtime_error("cufftCreate failed");

        std::unique_ptr<cufftHandle, std::function<void(cufftHandle*)>> const
                managed_plan(&unmanaged_plan, [](cufftHandle* handle) { cufftDestroy(*handle); });

        cufftSetStream(unmanaged_plan, stream);
        cufft_rt = cufftPlanMany(
                &unmanaged_plan, // plan handle
                sizeof...(X),
                n.data(), // Nx, Ny...
                NULL,
                1,
                idist,
                NULL,
                1,
                odist,
                cufft_transform_type<Tin, Tout>(),
                1);

        if (cufft_rt != CUFFT_SUCCESS)
            throw std::runtime_error("cufftPlan failed");

        cufft_rt = _cufftExec<Tin, Tout>(
                kwargs.direction == ddc::FFT_Direction::FORWARD ? CUFFT_FORWARD : CUFFT_INVERSE,
                unmanaged_plan,
                reinterpret_cast<typename _cufft_type<Tin>::type*>(in_data),
                reinterpret_cast<typename _cufft_type<Tout>::type*>(out_data));
        if (cufft_rt != CUFFT_SUCCESS)
            throw std::runtime_error("cufftExec failed");
        // std::cout << "performed with cufft";
    }
#endif
#if hipfft_AVAIL && !HIP_FOR_NVIDIA
    else if constexpr (std::is_same_v<ExecSpace, Kokkos::HIP>) {
        hipStream_t stream = execSpace.hip_stream();

        hipfftHandle unmanaged_plan;
        hipfftResult hipfft_rt = hipfftCreate(&unmanaged_plan);

        if (hipfft_rt != HIPFFT_SUCCESS)
            throw std::runtime_error("hipfftCreate failed");

        std::unique_ptr<hipfftHandle, std::function<void(hipfftHandle*)>> const
                managed_plan(&unmanaged_plan, [](hipfftHandle* handle) { hipfftDestroy(*handle); });

        hipfftSetStream(unmanaged_plan, stream);
        hipfft_rt = hipfftPlanMany(
                &unmanaged_plan, // plan handle
                sizeof...(X),
                n.data(), // Nx, Ny...
                NULL,
                1,
                idist,
                NULL,
                1,
                odist,
                hipfft_transform_type<Tin, Tout>(),
                1);

        if (hipfft_rt != HIPFFT_SUCCESS)
            throw std::runtime_error("hipfftPlan failed");

        hipfft_rt = _hipfftExec<Tin, Tout>(
                kwargs.direction == ddc::FFT_Direction::FORWARD ? HIPFFT_FORWARD : HIPFFT_BACKWARD,
                unmanaged_plan,
                reinterpret_cast<typename _hipfft_type<Tin>::type*>(in_data),
                reinterpret_cast<typename _hipfft_type<Tout>::type*>(out_data));
        if (hipfft_rt != HIPFFT_SUCCESS)
            throw std::runtime_error("hipfftExec failed");
        // std::cout << "performed with hipfft";
    }
#endif

    if (kwargs.normalization != ddc::FFT_Normalization::OFF) {
        real_type_t<Tout> norm_coef = 1;
        switch (kwargs.normalization) {
        case ddc::FFT_Normalization::OFF:
            break;
        case ddc::FFT_Normalization::FORWARD:
            norm_coef
                    = kwargs.direction == ddc::FFT_Direction::FORWARD
                              ? 1. / (ddc::get<ddc::UniformPointSampling<X>>(mesh.extents()) * ...)
                              : 1.;
            break;
        case ddc::FFT_Normalization::BACKWARD:
            norm_coef
                    = kwargs.direction == ddc::FFT_Direction::BACKWARD
                              ? 1. / (ddc::get<ddc::UniformPointSampling<X>>(mesh.extents()) * ...)
                              : 1.;
            break;
        case ddc::FFT_Normalization::ORTHO:
            norm_coef = 1.
                        / Kokkos::sqrt(
                                (ddc::get<ddc::UniformPointSampling<X>>(mesh.extents()) * ...));
            break;
        case ddc::FFT_Normalization::FULL:
            norm_coef
                    = kwargs.direction == ddc::FFT_Direction::FORWARD
                              ? (((coordinate(
                                           ddc::select<ddc::UniformPointSampling<X>>(mesh).back())
                                   - coordinate(
                                           ddc::select<ddc::UniformPointSampling<X>>(mesh).front()))
                                  / (ddc::get<ddc::UniformPointSampling<X>>(mesh.extents()) - 1)
                                  / Kokkos::sqrt(2 * Kokkos::numbers::pi))
                                 * ...)
                              : ((Kokkos::sqrt(2 * Kokkos::numbers::pi)
                                  / (coordinate(
                                             ddc::select<ddc::UniformPointSampling<X>>(mesh).back())
                                     - coordinate(ddc::select<ddc::UniformPointSampling<X>>(mesh)
                                                          .front()))
                                  * (ddc::get<ddc::UniformPointSampling<X>>(mesh.extents()) - 1)
                                  / ddc::get<ddc::UniformPointSampling<X>>(mesh.extents()))
                                 * ...);
            break;
        }

        Kokkos::parallel_for(
                Kokkos::RangePolicy<ExecSpace>(
                        execSpace,
                        0,
                        is_complex_v<Tout> && transform_type_v<Tin, Tout> != TransformType::C2C
                                ? (LastSelector<double, X, X...>(
                                           ddc::get<ddc::UniformPointSampling<X>>(mesh.extents())
                                                           / 2
                                                   + 1,
                                           ddc::get<ddc::UniformPointSampling<X>>(mesh.extents()))
                                   * ...)
                                : (ddc::get<ddc::UniformPointSampling<X>>(mesh.extents()) * ...)),
                KOKKOS_LAMBDA(const int& i) { out_data[i] = out_data[i] * norm_coef; });
    }
}

template <typename X>
typename ddc::PeriodicSampling<ddc::Fourier<X>>::template Impl<Kokkos::HostSpace> FourierSampling(
        ddc::DiscreteDomain<ddc::UniformPointSampling<X>> x_mesh)
{
    auto [impl, ddom] = ddc::PeriodicSampling<ddc::Fourier<X>>::
            init(ddc::Coordinate<ddc::Fourier<X>>(0),
                 ddc::Coordinate<ddc::Fourier<X>>(
                         2 * (ddc::detail::fft::N<X>(x_mesh) - 1)
                         / (ddc::detail::fft::b<X>(x_mesh) - ddc::detail::fft::a<X>(x_mesh))
                         * Kokkos::numbers::pi),
                 ddc::DiscreteVector<ddc::PeriodicSampling<ddc::Fourier<X>>>(
                         ddc::detail::fft::N<X>(x_mesh)),
                 ddc::DiscreteVector<ddc::PeriodicSampling<ddc::Fourier<X>>>(
                         ddc::detail::fft::N<X>(x_mesh)));
    return std::move(impl);
}
} // namespace ddc::detail::fft

namespace ddc {

template <typename... X>
void init_fourier_space(ddc::DiscreteDomain<ddc::UniformPointSampling<X>...> x_mesh)
{
    return (ddc::init_discrete_space<ddc::PeriodicSampling<ddc::Fourier<X>>>(
                    ddc::detail::fft::FourierSampling(
                            ddc::select<ddc::UniformPointSampling<X>>(x_mesh))),
            ...);
}

// FourierMesh, first element corresponds to mode 0
template <typename... X>
ddc::DiscreteDomain<ddc::PeriodicSampling<ddc::Fourier<X>>...> FourierMesh(
        ddc::DiscreteDomain<ddc::UniformPointSampling<X>...> x_mesh,
        bool C2C)
{
    return ddc::DiscreteDomain<ddc::PeriodicSampling<ddc::Fourier<X>>...>(
            ddc::DiscreteDomain<ddc::PeriodicSampling<ddc::Fourier<X>>>(
                    ddc::DiscreteElement<ddc::PeriodicSampling<ddc::Fourier<X>>>(0),
                    ddc::DiscreteVector<ddc::PeriodicSampling<ddc::Fourier<X>>>(
                            (C2C ? ddc::detail::fft::N<X>(x_mesh)
                                 : ddc::detail::fft::LastSelector<double, X, X...>(
                                         ddc::detail::fft::N<X>(x_mesh) / 2 + 1,
                                         ddc::detail::fft::N<X>(x_mesh)))))...);
}

struct kwArgs_fft
{
    ddc::FFT_Normalization normalization;
};

// FFT
template <
        typename Tin,
        typename Tout,
        typename... X,
        typename ExecSpace,
        typename MemorySpace,
        typename layout_in,
        typename layout_out>
void fft(
        ExecSpace execSpace,
        ddc::ChunkSpan<
                Tout,
                ddc::DiscreteDomain<ddc::PeriodicSampling<ddc::Fourier<X>>...>,
                layout_out,
                MemorySpace> out,
        ddc::ChunkSpan<
                Tin,
                ddc::DiscreteDomain<ddc::UniformPointSampling<X>...>,
                layout_in,
                MemorySpace> in,
        ddc::kwArgs_fft kwargs = {ddc::FFT_Normalization::OFF})
{
    static_assert(
            std::is_same_v<
                    layout_in,
                    std::experimental::
                            layout_right> && std::is_same_v<layout_out, std::experimental::layout_right>,
            "Layouts must be right-handed");

    ddc::detail::fft::core<Tin, Tout, ExecSpace, MemorySpace, X...>(
            execSpace,
            out.data_handle(),
            in.data_handle(),
            in.domain(),
            {ddc::FFT_Direction::FORWARD, kwargs.normalization});
}

// iFFT (deduced from the fact that "in" is identified as a function on the Fourier space)
template <
        typename Tin,
        typename Tout,
        typename... X,
        typename ExecSpace,
        typename MemorySpace,
        typename layout_in,
        typename layout_out>
void ifft(
        ExecSpace execSpace,
        ddc::ChunkSpan<
                Tout,
                ddc::DiscreteDomain<ddc::UniformPointSampling<X>...>,
                layout_out,
                MemorySpace> out,
        ddc::ChunkSpan<
                Tin,
                ddc::DiscreteDomain<ddc::PeriodicSampling<ddc::Fourier<X>>...>,
                layout_in,
                MemorySpace> in,
        ddc::kwArgs_fft kwargs = {ddc::FFT_Normalization::OFF})
{
    static_assert(
            std::is_same_v<
                    layout_in,
                    std::experimental::
                            layout_right> && std::is_same_v<layout_out, std::experimental::layout_right>,
            "Layouts must be right-handed");

    ddc::detail::fft::core<Tin, Tout, ExecSpace, MemorySpace, X...>(
            execSpace,
            out.data_handle(),
            in.data_handle(),
            out.domain(),
            {ddc::FFT_Direction::BACKWARD, kwargs.normalization});
}
} // namespace ddc
