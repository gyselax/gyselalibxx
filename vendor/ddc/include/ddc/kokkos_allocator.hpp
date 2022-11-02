// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <type_traits>

#include <Kokkos_Core.hpp>

namespace ddc {

template <class T, class MemorySpace>
class KokkosAllocator
{
public:
    using value_type = T;

    using memory_space = MemorySpace;

    template <class U>
    struct rebind
    {
        using other = KokkosAllocator<U, MemorySpace>;
    };

    constexpr KokkosAllocator() = default;

    constexpr KokkosAllocator(KokkosAllocator const& x) = default;

    constexpr KokkosAllocator(KokkosAllocator&& x) noexcept = default;

    template <class U>
    constexpr explicit KokkosAllocator(KokkosAllocator<U, MemorySpace> const&) noexcept
    {
    }

    ~KokkosAllocator() = default;

    constexpr KokkosAllocator& operator=(KokkosAllocator const& x) = default;

    constexpr KokkosAllocator& operator=(KokkosAllocator&& x) noexcept = default;

    template <class U>
    constexpr KokkosAllocator& operator=(KokkosAllocator<U, MemorySpace> const&) noexcept
    {
    }

    [[nodiscard]] T* allocate(std::size_t n) const
    {
        return static_cast<T*>(Kokkos::kokkos_malloc<MemorySpace>(sizeof(T) * n));
    }

    void deallocate(T* p, std::size_t) const
    {
        Kokkos::kokkos_free(p);
    }
};

template <class T, class MST, class U, class MSU>
constexpr bool operator==(KokkosAllocator<T, MST> const&, KokkosAllocator<U, MSU> const&) noexcept
{
    return std::is_same_v<KokkosAllocator<T, MST>, KokkosAllocator<U, MSU>>;
}

template <class T, class MST, class U, class MSU>
constexpr bool operator!=(KokkosAllocator<T, MST> const&, KokkosAllocator<U, MSU> const&) noexcept
{
    return !std::is_same_v<KokkosAllocator<T, MST>, KokkosAllocator<U, MSU>>;
}

template <class T>
using DeviceAllocator = KokkosAllocator<T, Kokkos::DefaultExecutionSpace::memory_space>;

template <class T>
using HostAllocator = KokkosAllocator<T, Kokkos::HostSpace>;

} // namespace ddc
