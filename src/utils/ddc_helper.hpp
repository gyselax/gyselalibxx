// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

namespace ddcHelper {
//TODO: this should be directly handled by ddc::Discretization really,
//      in the meantime, we do it ourselves
/**
 * Computes fluid moments of the distribution function 
 * Density, mean velocity and temperature.
 *
 * @param dom The domain on which the length should be calculated.
 *
 * @return The length of the domain.
 */
template <class IDim>
constexpr std::enable_if_t<!IDim::continuous_dimension_type::PERIODIC, double>
total_interval_length(ddc::DiscreteDomain<IDim> const& dom)
{
    return std::fabs(ddc::rlength(dom));
}

//TODO: this should be directly handled by ddc::Discretization really,
//      in the meantime, we do it ourselves
/**
 * Computes fluid moments of the distribution function 
 * Density, mean velocity and temperature.
 *
 * @param dom The domain on which the length should be calculated.
 *
 * @return The length of the domain.
 */
template <class RDim>
constexpr std::enable_if_t<RDim::PERIODIC, double> total_interval_length(
        ddc::DiscreteDomain<ddc::UniformPointSampling<RDim>> const& dom)
{
    return std::fabs(ddc::rlength(dom) + ddc::step<ddc::UniformPointSampling<RDim>>());
}
}; // namespace ddcHelper

namespace detail {

/// \cond
template <class>
struct Device
{
};

template <class ElementType, class SupportType, class Allocator>
struct Device<ddc::Chunk<ElementType, SupportType, Allocator>>
{
    using type = typename ddc::Chunk<
            ElementType,
            SupportType,
            ddc::KokkosAllocator<ElementType, Kokkos::DefaultExecutionSpace::memory_space>>;
};

template <class ElementType, class SupportType, class Layout, class MemorySpace>
struct Device<ddc::ChunkSpan<ElementType, SupportType, Layout, MemorySpace>>
{
    using type = typename ddc::ChunkSpan<
            ElementType,
            SupportType,
            Layout,
            Kokkos::DefaultExecutionSpace::memory_space>;
};
/// \endcond

} // namespace detail

/// Alias template helper returning the "device" version of a `ddc::Chunk` or a `ddc::Chunkspan`
template <class C>
using device_t = typename detail::Device<C>::type;
