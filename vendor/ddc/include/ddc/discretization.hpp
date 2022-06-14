// SPDX-License-Identifier: MIT

#pragma once

#include <stdexcept>

#include <Kokkos_Core.hpp>
#if defined(__CUDACC__)
#include <cuda.h>
#endif
#if defined(__HIPCC__)
#include <hip/hip_runtime.h>
#endif

#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_domain.hpp"

namespace detail {

// For now, in the future, this should be specialized by tag
template <class IDimImpl>
inline IDimImpl* discretization_host = nullptr;

#if defined(__CUDACC__) || defined(__HIPCC__)
template <class IDimImpl>
inline __device__ __constant__ IDimImpl* discretization_device = nullptr;
#endif

template <class Tuple, std::size_t... Ids>
auto extract_after(Tuple&& t, std::index_sequence<Ids...>)
{
    return std::make_tuple(std::move(std::get<Ids + 1>(t))...);
}

template <class IDim>
void init_discretization_devices()
{
    using IDimImplHost = typename IDim::template Impl<Kokkos::HostSpace>;
#if defined(__CUDACC__)
    using IDimImplDevice = typename IDim::template Impl<Kokkos::CudaSpace>;
    discretization_host<IDimImplDevice> = new IDimImplDevice(*discretization_host<IDimImplHost>);
    IDimImplDevice* ptr_device;
    cudaMalloc(&ptr_device, sizeof(IDimImplDevice));
    cudaMemcpy(
            (void*)ptr_device,
            discretization_host<IDimImplDevice>,
            sizeof(IDimImplDevice),
            cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(
            discretization_device<IDimImplDevice>,
            &ptr_device,
            sizeof(IDimImplDevice*),
            0,
            cudaMemcpyHostToDevice);
#endif
#if defined(__HIPCC__)
    using IDimImplDevice = typename IDim::template Impl<Kokkos::Experimental::HIPSpace>;
    discretization_host<IDimImplDevice> = new IDimImplDevice(*discretization_host<IDimImplHost>);
    IDimImplDevice* ptr_device;
    hipMalloc(&ptr_device, sizeof(IDimImplDevice));
    hipMemcpy(
            (void*)ptr_device,
            discretization_host<IDimImplDevice>,
            sizeof(IDimImplDevice),
            hipMemcpyHostToDevice);
    hipMemcpyToSymbol(
            discretization_device<IDimImplDevice>,
            &ptr_device,
            sizeof(IDimImplDevice*),
            0,
            hipMemcpyHostToDevice);
#endif
}

} // namespace detail

template <class IDimImpl, class Arg>
Arg init_discretization(std::tuple<IDimImpl, Arg>&& a)
{
    using IDim = typename IDimImpl::ddim_type;
    using IDimImplHost = typename IDim::template Impl<Kokkos::HostSpace>;
    if (detail::discretization_host<IDimImplHost>) {
        throw std::runtime_error("Discretization function already initialized.");
    }
    detail::discretization_host<IDimImplHost> = new IDimImplHost(std::move(std::get<0>(a)));
    detail::init_discretization_devices<IDim>();
    return std::get<1>(a);
}

template <class IDimImpl, class... Args>
std::enable_if_t<2 <= sizeof...(Args), std::tuple<Args...>> init_discretization(
        std::tuple<IDimImpl, Args...>&& a)
{
    using IDim = typename IDimImpl::ddim_type;
    using IDimImplHost = typename IDim::template Impl<Kokkos::HostSpace>;
    if (detail::discretization_host<IDimImplHost>) {
        throw std::runtime_error("Discretization function already initialized.");
    }
    detail::discretization_host<IDimImplHost> = new IDimImplHost(std::move(std::get<0>(a)));
    detail::init_discretization_devices<IDim>();
    return detail::extract_after(std::move(a), std::index_sequence_for<Args...>());
}

template <class IDim, class... Args>
void init_discretization(Args&&... a)
{
    using IDimImplHost = typename IDim::template Impl<Kokkos::HostSpace>;
    if (detail::discretization_host<std::remove_cv_t<std::remove_reference_t<IDimImplHost>>>) {
        throw std::runtime_error("Discretization function already initialized.");
    }
    detail::discretization_host<IDimImplHost> = new IDimImplHost(std::forward<Args>(a)...);
    detail::init_discretization_devices<IDim>();
}

template <class IDim>
inline typename IDim::template Impl<Kokkos::HostSpace> const& discretization()
{
    return *detail::discretization_host<typename IDim::template Impl<Kokkos::HostSpace>>;
}

template <class IDim>
inline typename IDim::template Impl<Kokkos::HostSpace> const& discretization_host()
{
    return *detail::discretization_host<typename IDim::template Impl<Kokkos::HostSpace>>;
}

#if defined(__CUDACC__)
template <class IDim>
__device__ inline typename IDim::template Impl<Kokkos::CudaSpace> const& discretization_device()
{
    return *detail::discretization_device<typename IDim::template Impl<Kokkos::CudaSpace>>;
}
#endif

#if defined(__HIPCC__)
template <class IDim>
__device__ inline typename IDim::template Impl<Kokkos::Experimental::HIPSpace> const& discretization_device()
{
    return *detail::discretization_device<typename IDim::template Impl<Kokkos::Experimental::HIPSpace>>;
}
#endif
