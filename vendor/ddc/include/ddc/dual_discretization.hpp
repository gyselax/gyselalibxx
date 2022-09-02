// SPDX-License-Identifier: MIT

#pragma once

#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

#include <Kokkos_Core.hpp>

#include "ddc/detail/macros.hpp"

#if defined(__CUDACC__)
#include <cuda.h>
#endif
#if defined(__HIPCC__)
#include <hip/hip_runtime.h>
#endif

template <class DDim>
class DualDiscretization
{
    using DDimImplHost = typename DDim::template Impl<Kokkos::HostSpace>;
#if defined(__CUDACC__)
    using DDimImplDevice = typename DDim::template Impl<Kokkos::CudaSpace>;
#elif defined(__HIPCC__)
    using DDimImplDevice = typename DDim::template Impl<Kokkos::Experimental::HIPSpace>;
#endif

    DDimImplHost m_host;
#if defined(__CUDACC__) || defined(__HIPCC__)
    DDimImplDevice m_device_on_host;
    std::unique_ptr<DDimImplDevice, std::function<void(DDimImplDevice*)>> m_device;
#endif

public:
    template <class... Args>
    explicit DualDiscretization(Args&&... args)
        : m_host(std::forward<Args>(args)...)
#if defined(__CUDACC__) || defined(__HIPCC__)
        , m_device_on_host(m_host)
#endif
    {
#if defined(__CUDACC__)
        DDimImplDevice* ptr_device;
        cudaMalloc(&ptr_device, sizeof(DDimImplDevice));
        m_device = std::unique_ptr<
                DDimImplDevice,
                std::function<void(DDimImplDevice*)>>(ptr_device, [](DDimImplDevice* ptr) {
            cudaFree(ptr);
        });
        cudaMemcpy(
                reinterpret_cast<void*>(ptr_device),
                &m_device_on_host,
                sizeof(DDimImplDevice),
                cudaMemcpyHostToDevice);
#elif defined(__HIPCC__)
        DDimImplDevice* ptr_device;
        hipMalloc(&ptr_device, sizeof(DDimImplDevice));
        m_device = std::unique_ptr<
                DDimImplDevice,
                std::function<void(DDimImplDevice*)>>(ptr_device, [](DDimImplDevice* ptr) {
            hipFree(ptr);
        });
        hipMemcpy(
                reinterpret_cast<void*>(ptr_device),
                &m_device_on_host,
                sizeof(DDimImplDevice),
                hipMemcpyHostToDevice);
#endif
    }

    template <class MemorySpace>
    DDC_INLINE_FUNCTION typename DDim::template Impl<MemorySpace> const& get()
    {
        if constexpr (std::is_same_v<MemorySpace, Kokkos::HostSpace>) {
            return m_host;
        }
#if defined(__CUDACC__)
        else if constexpr (std::is_same_v<MemorySpace, Kokkos::CudaSpace>) {
            return m_device_on_host;
        }
#elif defined(__CUDACC__)
        else if constexpr (std::is_same_v<MemorySpace, Kokkos::Experimental::HIPSpace>) {
            return m_device_on_host;
        }
#endif
        else {
            static_assert(!std::is_same_v<MemorySpace, MemorySpace>);
        }
    }

#if defined(__CUDACC__) || defined(__HIPCC__)
    DDimImplDevice* get_device_ptr() const
    {
        return m_device.get();
    }
#endif
};
