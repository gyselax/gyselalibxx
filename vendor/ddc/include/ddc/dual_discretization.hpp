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

namespace ddc {

template <class DDim>
class DualDiscretization
{
    using DDimImplHost = typename DDim::template Impl<Kokkos::HostSpace>;
#if defined(__CUDACC__)
    using DDimImplDevice = typename DDim::template Impl<Kokkos::CudaSpace>;
#elif defined(__HIPCC__)
    using DDimImplDevice = typename DDim::template Impl<Kokkos::HIPSpace>;
#else
    using DDimImplDevice = DDimImplHost;
#endif

    DDimImplHost m_host;
#if defined(__CUDACC__) || defined(__HIPCC__)
    DDimImplDevice m_device_on_host;
#endif

public:
    template <class... Args>
    explicit DualDiscretization(Args&&... args)
        : m_host(std::forward<Args>(args)...)
#if defined(__CUDACC__) || defined(__HIPCC__)
        , m_device_on_host(m_host)
#endif
    {
    }

    template <class MemorySpace>
    KOKKOS_FORCEINLINE_FUNCTION typename DDim::template Impl<MemorySpace> const& get()
    {
        if constexpr (std::is_same_v<MemorySpace, Kokkos::HostSpace>) {
            return m_host;
        }
#if defined(__CUDACC__)
        else if constexpr (std::is_same_v<MemorySpace, Kokkos::CudaSpace>) {
            return m_device_on_host;
        }
#elif defined(__HIPCC__)
        else if constexpr (std::is_same_v<MemorySpace, Kokkos::HIPSpace>) {
            return m_device_on_host;
        }
#endif
        else {
            static_assert(!std::is_same_v<MemorySpace, MemorySpace>);
        }
    }

    KOKKOS_FORCEINLINE_FUNCTION DDimImplHost const& get_host()
    {
        return m_host;
    }

    KOKKOS_FORCEINLINE_FUNCTION DDimImplDevice const& get_device()
    {
#if defined(__CUDACC__) || defined(__HIPCC__)
        return m_device_on_host;
#else
        return m_host;
#endif
    }
};

} // namespace ddc
