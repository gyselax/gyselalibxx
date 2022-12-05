// SPDX-License-Identifier: MIT

#pragma once

#include <map>
#include <memory>
#include <optional>
#include <ostream>
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
#include "ddc/dual_discretization.hpp"

namespace ddc {

namespace ddc_detail {

template <class DDim, class MemorySpace>
using ddim_impl_t = typename DDim::template Impl<MemorySpace>;

template <class T>
class gpu_proxy
{
    // Here are some reasonable concepts that T should satisfy to avoid undefined behaviors:
    // - copy-constructible: objects may be memcopied to the device,
    // - standard layout: objects will be ensured to have the same, standard, representation on the host and the device,
    // - trivially destructible: the destructor of objects located on a device may not be called.
    // static_assert(std::is_standard_layout_v<T>, "Not standard layout");
    // static_assert(std::is_trivially_destructible_v<T>, "Not trivially destructible");
    // static_assert(std::is_trivially_copy_constructible_v<T>, "Not trivially copy-constructible");
    // Currently not trivially destructible because for example of the Kokkos::View (mostly reference-counting)
    // Currently not trivially copy-constructible because of discrete spaces that have deleted copy-constructors and Kokkos::View (mostly reference-counting)

private:
    alignas(T) std::byte m_data[sizeof(T)];

public:
    DDC_INLINE_FUNCTION
    T* operator->()
    {
        return reinterpret_cast<T*>(m_data);
    }

    DDC_INLINE_FUNCTION
    T& operator*()
    {
        return *reinterpret_cast<T*>(m_data);
    }

    DDC_INLINE_FUNCTION
    T* data()
    {
        return reinterpret_cast<T*>(m_data);
    }
};

// Global CPU variable storing resetters. Required to correctly free data.
inline std::optional<std::map<std::string, std::function<void()>>> g_discretization_store;

// Global CPU variable owning discrete spaces data for CPU and GPU
template <class DDim>
inline std::optional<DualDiscretization<DDim>> g_discrete_space_dual;

#if defined(__CUDACC__)
// Global GPU variable viewing data owned by the CPU
template <class DDim>
__constant__ gpu_proxy<ddim_impl_t<DDim, Kokkos::CudaSpace>> g_discrete_space_device;
#elif defined(__HIPCC__)
// Global GPU variable viewing data owned by the CPU
// WARNING: do not put the `inline` keyword, seems to fail on MI100 rocm/4.5.0
template <class DDim>
__constant__ gpu_proxy<ddim_impl_t<DDim, Kokkos::Experimental::HIPSpace>> g_discrete_space_device;
#endif

inline void display_discretization_store(std::ostream& os)
{
    if (g_discretization_store) {
        os << "The host discretization store is initialized:\n";
        for (auto const& [key, value] : *g_discretization_store) {
            os << " - " << key << "\n";
        }
    } else {
        os << "The host discretization store is not initialized:\n";
    }
}

template <class Tuple, std::size_t... Ids>
auto extract_after(Tuple&& t, std::index_sequence<Ids...>)
{
    return std::make_tuple(std::move(std::get<Ids + 1>(t))...);
}

} // namespace ddc_detail

/** Initialize (emplace) a global singleton discrete space
 * 
 * @param args the constructor arguments
 */
template <class DDim, class... Args>
void init_discrete_space(Args&&... args)
{
    if (ddc_detail::g_discrete_space_dual<DDim>) {
        throw std::runtime_error("Discrete space function already initialized.");
    }
    ddc_detail::g_discrete_space_dual<DDim>.emplace(std::forward<Args>(args)...);
    ddc_detail::g_discretization_store->emplace(typeid(DDim).name(), []() {
        ddc_detail::g_discrete_space_dual<DDim>.reset();
    });
#if defined(__CUDACC__)
    cudaMemcpyToSymbol(
            ddc_detail::g_discrete_space_device<DDim>,
            &ddc_detail::g_discrete_space_dual<DDim>->get_device(),
            sizeof(ddc_detail::g_discrete_space_dual<DDim>->get_device()));
#elif defined(__HIPCC__)
    hipMemcpyToSymbol(
            ddc_detail::g_discrete_space_device<DDim>,
            &ddc_detail::g_discrete_space_dual<DDim>->get_device(),
            sizeof(ddc_detail::g_discrete_space_dual<DDim>->get_device()));
#endif
}

/** Move construct a global singleton discrete space and pass through the other argument
 * 
 * @param a - the discrete space to move at index 0
 *          - the arguments to pass through at index 1
 */
template <class DDimImpl, class Arg>
Arg init_discrete_space(std::tuple<DDimImpl, Arg>&& a)
{
    using DDim = typename DDimImpl::discrete_dimension_type;
    init_discrete_space<DDim>(std::move(std::get<0>(a)));
    return std::get<1>(a);
}

/** Move construct a global singleton discrete space and pass through remaining arguments
 * 
 * @param a - the discrete space to move at index 0
 *          - the (2+) arguments to pass through in other indices
 */
template <class DDimImpl, class... Args>
std::enable_if_t<2 <= sizeof...(Args), std::tuple<Args...>> init_discrete_space(
        std::tuple<DDimImpl, Args...>&& a)
{
    using DDim = typename DDimImpl::discrete_dimension_type;
    init_discrete_space<DDim>(std::move(std::get<0>(a)));
    return ddc_detail::extract_after(std::move(a), std::index_sequence_for<Args...>());
}

template <class DDim, class MemorySpace = DDC_CURRENT_KOKKOS_SPACE>
DDC_INLINE_FUNCTION ddc_detail::ddim_impl_t<DDim, MemorySpace> const& discrete_space()
{
    if constexpr (std::is_same_v<MemorySpace, Kokkos::HostSpace>) {
        return ddc_detail::g_discrete_space_dual<DDim>->get_host();
    }
#if defined(__CUDACC__)
    else if constexpr (std::is_same_v<MemorySpace, Kokkos::CudaSpace>) {
        return *ddc_detail::g_discrete_space_device<DDim>;
    }
#elif defined(__HIPCC__)
    else if constexpr (std::is_same_v<MemorySpace, Kokkos::Experimental::HIPSpace>) {
        return *ddc_detail::g_discrete_space_device<DDim>;
    }
#endif
    else {
        static_assert(std::is_same_v<MemorySpace, MemorySpace>, "Memory space not handled");
    }
}

} // namespace ddc
