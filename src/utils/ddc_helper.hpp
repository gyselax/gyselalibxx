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
