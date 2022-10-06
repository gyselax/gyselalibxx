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

namespace detail {

template <class DDim, class MemorySpace>
struct DiscreteSpaceGetter;

inline std::optional<std::map<std::string, std::function<void()>>> g_discretization_store;

// For now, in the future, this should be specialized by tag
template <class DDim>
inline std::optional<DualDiscretization<DDim>> g_discrete_space_host;

template <class DDim>
struct DiscreteSpaceGetter<DDim, Kokkos::HostSpace>
{
    static inline typename DDim::template Impl<Kokkos::HostSpace> const& get()
    {
        return g_discrete_space_host<DDim>->template get<Kokkos::HostSpace>();
    }
};

#if defined(__CUDACC__) || defined(__HIPCC__)
// WARNING: do not put the `inline` keyword, seems to fail on MI100 rocm/4.5.0
template <class DDimImpl>
__device__ __constant__ DDimImpl* g_discrete_space_device = nullptr;

template <class DDim, class MemorySpace>
struct DiscreteSpaceGetter
{
    DDC_INLINE_FUNCTION
    static typename DDim::template Impl<MemorySpace> const& get()
    {
        return *g_discrete_space_device<typename DDim::template Impl<MemorySpace>>;
    }
};
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

} // namespace detail

/** Initialize (emplace) a global singleton discrete space
 * 
 * @param a the constructor arguments
 */
template <class DDim, class... Args>
void init_discrete_space(Args&&... args)
{
    if (detail::g_discrete_space_host<DDim>) {
        throw std::runtime_error("Discrete space function already initialized.");
    }
    detail::g_discrete_space_host<DDim>.emplace(std::forward<Args>(args)...);
    detail::g_discretization_store->emplace(typeid(DDim).name(), []() {
        detail::g_discrete_space_host<DDim>.reset();
    });
#if defined(__CUDACC__)
    using DDimImplDevice = typename DDim::template Impl<Kokkos::CudaSpace>;
    DDimImplDevice* device_ptr = detail::g_discrete_space_host<DDim>->get_device_ptr();
    cudaMemcpyToSymbol(
            detail::g_discrete_space_device<DDimImplDevice>,
            &device_ptr,
            sizeof(DDimImplDevice*),
            0,
            cudaMemcpyHostToDevice);
#endif
#if defined(__HIPCC__)
    using DDimImplDevice = typename DDim::template Impl<Kokkos::Experimental::HIPSpace>;
    DDimImplDevice* device_ptr = detail::g_discrete_space_host<DDim>->get_device_ptr();
    hipMemcpyToSymbol(
            detail::g_discrete_space_device<DDimImplDevice>,
            &device_ptr,
            sizeof(DDimImplDevice*),
            0,
            hipMemcpyHostToDevice);
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
    return detail::extract_after(std::move(a), std::index_sequence_for<Args...>());
}

template <class DDim, class MemorySpace = DDC_CURRENT_KOKKOS_SPACE>
DDC_INLINE_FUNCTION typename DDim::template Impl<MemorySpace> const& discrete_space()
{
    return detail::DiscreteSpaceGetter<DDim, MemorySpace>::get();
}
