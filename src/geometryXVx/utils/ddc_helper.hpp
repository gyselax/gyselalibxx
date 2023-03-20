// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

/**
 * Computes fluid moments of the distribution function 
 * Density, mean velocity and temperature.
 */
class ddcHelper
{
public:
    //TODO: this should be directly handled by ddc::Discretization really,
    //      in the meantime, we do it ourselves
    template <class IDim>
    static constexpr std::enable_if_t<!IDim::continuous_dimension_type::PERIODIC, double>
    total_interval_length(ddc::DiscreteDomain<IDim> const& dom)
    {
        return std::fabs(ddc::rlength(dom));
    }

    //TODO: this should be directly handled by ddc::Discretization really,
    //      in the meantime, we do it ourselves
    template <class RDim>
    static constexpr std::enable_if_t<RDim::PERIODIC, double> total_interval_length(
            ddc::DiscreteDomain<ddc::UniformPointSampling<RDim>> const& dom)
    {
        return std::fabs(ddc::rlength(dom) + ddc::step<ddc::UniformPointSampling<RDim>>());
    }
};
