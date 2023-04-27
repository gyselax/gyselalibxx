// SPDX-License-Identifier: MIT

#pragma once

#include "ddc/coordinate.hpp"
#include "ddc/discrete_element.hpp"
#include "ddc/discrete_space.hpp"
#include "ddc/non_uniform_point_sampling.hpp"
#include "ddc/uniform_point_sampling.hpp"

namespace ddc {

template <
        class... PointSampling,
        std::enable_if_t<
                ((is_uniform_sampling_v<
                          PointSampling> || is_non_uniform_sampling_v<PointSampling>)&&...),
                int> = 0>
DDC_INLINE_FUNCTION Coordinate<typename PointSampling::continuous_dimension_type...> coordinate(
        DiscreteElement<PointSampling...> const& c)
{
    return Coordinate<typename PointSampling::continuous_dimension_type...>(
            discrete_space<PointSampling>().coordinate(select<PointSampling>(c))...);
}

} // namespace ddc
