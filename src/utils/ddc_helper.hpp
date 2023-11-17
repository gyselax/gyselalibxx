// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

namespace ddcHelper {
//TODO: this should be directly handled by ddc::Discretization really,
//      in the meantime, we do it ourselves
/**
 * Calculate the total length of a non-periodic domain.
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
 * Calculate the total length of a uniform periodic domain.
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

//TODO: this should be directly handled by ddc::Discretization really,
//      in the meantime, we do it ourselves
/**
 * Calculate the total length of a non-uniform periodic domain.
 *
 * @param dom The domain on which the length should be calculated.
 *
 * @return The length of the domain.
 */
template <class RDim>
constexpr std::enable_if_t<RDim::PERIODIC, double> total_interval_length(
        ddc::DiscreteDomain<ddc::NonUniformPointSampling<RDim>> const& dom)
{
    using IDim = ddc::NonUniformPointSampling<RDim>;
    ddc::DiscreteDomain<IDim> dom_periodic(dom.front(), dom.extents() + 1);
    return std::fabs(ddc::rlength(dom_periodic));
}

//TODO: this should be directly handled by ddc::Discretization really,
//      in the meantime, we do it ourselves
/**
 * @brief Calculate the Coordinate inside the domain.
 *
 * In the case of a periodic domain, a Coordinate can sometimes be found
 * outside the domain. In this case it is useful to be able to find the
 * equivalent coordinate inside the domain. This function makes that
 * possible.
 *
 * @param[in] coord
 *      The 1D coordinate we want to compute inside the domain.
 * @param[in] dom
 *      The domain where the coordinate is defined.
 *
 * @return The equivalent coordinate inside the domain.
 */
template <class IDim>
constexpr std::enable_if_t<
        IDim::continuous_dimension_type::PERIODIC,
        typename IDim::continuous_element_type>
restrict_to_domain(
        typename IDim::continuous_element_type coord,
        ddc::DiscreteDomain<IDim> const& dom)
{
    using Coord = typename IDim::continuous_element_type;
    double const x_min = ddc::rmin(dom);
    double const length = total_interval_length(dom);
    double const x_max = x_min + length;

    assert(length > 0);
    coord -= x_min;
    if (fabs(coord) > 10 * length) {
        double periodic_factor = 2 * M_PI / length;
        double coord_2pi = double(coord) * periodic_factor;
        coord_2pi = std::copysign(std::acos(std::cos(coord_2pi)), std::sin(coord_2pi));
        coord = coord_2pi < 0 ? Coord((coord_2pi + 2 * M_PI) / periodic_factor)
                              : Coord((coord_2pi) / periodic_factor);
    }
    coord += x_min;
    while (coord < x_min)
        coord += length;
    while (coord >= x_max)
        coord -= length;
    return coord;
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
