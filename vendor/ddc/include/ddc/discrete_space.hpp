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

template <class IDim, class MemorySpace>
struct DiscreteSpaceGetter;

// For now, in the future, this should be specialized by tag
template <class IDimImpl>
inline IDimImpl* g_discrete_space_host = nullptr;

template <class IDim>
struct DiscreteSpaceGetter<IDim, Kokkos::HostSpace>
{
    static inline typename IDim::template Impl<Kokkos::HostSpace> const& get()
    {
        return *g_discrete_space_host<typename IDim::template Impl<Kokkos::HostSpace>>;
    }
};

#if defined(__CUDACC__) || defined(__HIPCC__)
// WARNING: do not put the `inline` keyword, seems to fail on MI100 rocm/4.5.0
template <class IDimImpl>
__device__ __constant__ IDimImpl* g_discrete_space_device = nullptr;

template <class IDim, class MemorySpace>
struct DiscreteSpaceGetter
{
    DDC_INLINE_FUNCTION
    static typename IDim::template Impl<MemorySpace> const& get()
    {
        return *g_discrete_space_device<typename IDim::template Impl<MemorySpace>>;
    }
};
#endif

template <class Tuple, std::size_t... Ids>
auto extract_after(Tuple&& t, std::index_sequence<Ids...>)
{
    return std::make_tuple(std::move(std::get<Ids + 1>(t))...);
}

template <class IDim>
void init_discrete_space_devices()
{
#if defined(__CUDACC__)
    using IDimImplHost = typename IDim::template Impl<Kokkos::HostSpace>;
    using IDimImplDevice = typename IDim::template Impl<Kokkos::CudaSpace>;
    g_discrete_space_host<IDimImplDevice> = new IDimImplDevice(
            *g_discrete_space_host<IDimImplHost>);
    IDimImplDevice* ptr_device;
    cudaMalloc(&ptr_device, sizeof(IDimImplDevice));
    cudaMemcpy(
            (void*)ptr_device,
            g_discrete_space_host<IDimImplDevice>,
            sizeof(IDimImplDevice),
            cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(
            g_discrete_space_device<IDimImplDevice>,
            &ptr_device,
            sizeof(IDimImplDevice*),
            0,
            cudaMemcpyHostToDevice);
#endif
#if defined(__HIPCC__)
    using IDimImplHost = typename IDim::template Impl<Kokkos::HostSpace>;
    using IDimImplDevice = typename IDim::template Impl<Kokkos::Experimental::HIPSpace>;
    g_discrete_space_host<IDimImplDevice> = new IDimImplDevice(
            *g_discrete_space_host<IDimImplHost>);
    IDimImplDevice* ptr_device;
    hipMalloc(&ptr_device, sizeof(IDimImplDevice));
    hipMemcpy(
            (void*)ptr_device,
            g_discrete_space_host<IDimImplDevice>,
            sizeof(IDimImplDevice),
            hipMemcpyHostToDevice);
    hipMemcpyToSymbol(
            g_discrete_space_device<IDimImplDevice>,
            &ptr_device,
            sizeof(IDimImplDevice*),
            0,
            hipMemcpyHostToDevice);
#endif
}

} // namespace detail

/** Initialize (emplace) a global singleton discrete space
 * 
 * @param a the constructor arguments
 */
template <class DDim, class... Args>
static inline void init_discrete_space(Args&&... a)
{
    using DDimImplHost = typename DDim::template Impl<Kokkos::HostSpace>;
    if (detail::g_discrete_space_host<std::remove_cv_t<std::remove_reference_t<DDimImplHost>>>) {
        throw std::runtime_error("Discrete space function already initialized.");
    }
    detail::g_discrete_space_host<DDimImplHost> = new DDimImplHost(std::forward<Args>(a)...);
    detail::init_discrete_space_devices<DDim>();
}

/** Move construct a global singleton discrete space and pass through the other argument
 * 
 * @param a - the discrete space to move at index 0
 *          - the arguments to pass through at index 1
 */
template <class IDimImpl, class Arg>
static inline Arg init_discrete_space(std::tuple<IDimImpl, Arg>&& a)
{
    using DDim = typename IDimImpl::discrete_dimension_type;
    init_discrete_space<DDim>(std::move(std::get<0>(a)));
    return std::get<1>(a);
}

/** Move construct a global singleton discrete space and pass through remaining arguments
 * 
 * @param a - the discrete space to move at index 0
 *          - the (2+) arguments to pass through in other indices
 */
template <class IDimImpl, class... Args>
static inline std::enable_if_t<2 <= sizeof...(Args), std::tuple<Args...>> init_discrete_space(
        std::tuple<IDimImpl, Args...>&& a)
{
    using DDim = typename IDimImpl::discrete_dimension_type;
    init_discrete_space<DDim>(std::move(std::get<0>(a)));
    return detail::extract_after(std::move(a), std::index_sequence_for<Args...>());
}

template <class IDim, class MemorySpace = DDC_CURRENT_KOKKOS_SPACE>
DDC_INLINE_FUNCTION typename IDim::template Impl<MemorySpace> const& discrete_space()
{
    return detail::DiscreteSpaceGetter<IDim, MemorySpace>::get();
}
