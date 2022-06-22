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

#include "ddc/discrete_domain.hpp"
#include "ddc/discrete_space.hpp"

namespace detail {

// For now, in the future, this should be specialized by tag
template <class IDimImpl>
inline IDimImpl* discrete_space_host = nullptr;

#if defined(__CUDACC__) || defined(__HIPCC__)
// WARNING: do not put the `inline` keyword, seems to fail on MI100 rocm/4.5.0
template <class IDimImpl>
__device__ __constant__ IDimImpl* discrete_space_device = nullptr;
#endif

template <class Tuple, std::size_t... Ids>
auto extract_after(Tuple&& t, std::index_sequence<Ids...>)
{
    return std::make_tuple(std::move(std::get<Ids + 1>(t))...);
}

template <class IDim>
void init_discrete_space_devices()
{
    using IDimImplHost = typename IDim::template Impl<Kokkos::HostSpace>;
#if defined(__CUDACC__)
    using IDimImplDevice = typename IDim::template Impl<Kokkos::CudaSpace>;
    discrete_space_host<IDimImplDevice> = new IDimImplDevice(*discrete_space_host<IDimImplHost>);
    IDimImplDevice* ptr_device;
    cudaMalloc(&ptr_device, sizeof(IDimImplDevice));
    cudaMemcpy(
            (void*)ptr_device,
            discrete_space_host<IDimImplDevice>,
            sizeof(IDimImplDevice),
            cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(
            discrete_space_device<IDimImplDevice>,
            &ptr_device,
            sizeof(IDimImplDevice*),
            0,
            cudaMemcpyHostToDevice);
#endif
#if defined(__HIPCC__)
    using IDimImplDevice = typename IDim::template Impl<Kokkos::Experimental::HIPSpace>;
    discrete_space_host<IDimImplDevice> = new IDimImplDevice(*discrete_space_host<IDimImplHost>);
    IDimImplDevice* ptr_device;
    hipMalloc(&ptr_device, sizeof(IDimImplDevice));
    hipMemcpy(
            (void*)ptr_device,
            discrete_space_host<IDimImplDevice>,
            sizeof(IDimImplDevice),
            hipMemcpyHostToDevice);
    hipMemcpyToSymbol(
            discrete_space_device<IDimImplDevice>,
            &ptr_device,
            sizeof(IDimImplDevice*),
            0,
            hipMemcpyHostToDevice);
#endif
}

} // namespace detail

template <class IDimImpl, class Arg>
Arg init_discrete_space(std::tuple<IDimImpl, Arg>&& a)
{
    using IDim = typename IDimImpl::discrete_dimension_type;
    using IDimImplHost = typename IDim::template Impl<Kokkos::HostSpace>;
    if (detail::discrete_space_host<IDimImplHost>) {
        throw std::runtime_error("Discrete space function already initialized.");
    }
    detail::discrete_space_host<IDimImplHost> = new IDimImplHost(std::move(std::get<0>(a)));
    detail::init_discrete_space_devices<IDim>();
    return std::get<1>(a);
}

template <class IDimImpl, class... Args>
std::enable_if_t<2 <= sizeof...(Args), std::tuple<Args...>> init_discrete_space(
        std::tuple<IDimImpl, Args...>&& a)
{
    using IDim = typename IDimImpl::discrete_dimension_type;
    using IDimImplHost = typename IDim::template Impl<Kokkos::HostSpace>;
    if (detail::discrete_space_host<IDimImplHost>) {
        throw std::runtime_error("Discrete space function already initialized.");
    }
    detail::discrete_space_host<IDimImplHost> = new IDimImplHost(std::move(std::get<0>(a)));
    detail::init_discrete_space_devices<IDim>();
    return detail::extract_after(std::move(a), std::index_sequence_for<Args...>());
}

template <class IDim, class... Args>
void init_discrete_space(Args&&... a)
{
    using IDimImplHost = typename IDim::template Impl<Kokkos::HostSpace>;
    if (detail::discrete_space_host<std::remove_cv_t<std::remove_reference_t<IDimImplHost>>>) {
        throw std::runtime_error("Discrete space function already initialized.");
    }
    detail::discrete_space_host<IDimImplHost> = new IDimImplHost(std::forward<Args>(a)...);
    detail::init_discrete_space_devices<IDim>();
}

template <class IDim>
inline typename IDim::template Impl<Kokkos::HostSpace> const& discrete_space()
{
    return *detail::discrete_space_host<typename IDim::template Impl<Kokkos::HostSpace>>;
}

template <class IDim>
inline typename IDim::template Impl<Kokkos::HostSpace> const& discrete_space_host()
{
    return *detail::discrete_space_host<typename IDim::template Impl<Kokkos::HostSpace>>;
}

#if defined(__CUDACC__)
template <class IDim>
__device__ inline typename IDim::template Impl<Kokkos::CudaSpace> const& discrete_space_device()
{
    return *detail::discrete_space_device<typename IDim::template Impl<Kokkos::CudaSpace>>;
}
#endif

#if defined(__HIPCC__)
template <class IDim>
__device__ inline typename IDim::template Impl<Kokkos::Experimental::HIPSpace> const&
discrete_space_device()
{
    return *detail::discrete_space_device<
            typename IDim::template Impl<Kokkos::Experimental::HIPSpace>>;
}
#endif
