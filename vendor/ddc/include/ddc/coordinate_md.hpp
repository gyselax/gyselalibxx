// SPDX-License-Identifier: MIT

#pragma once

#include "ddc/coordinate.hpp"
#include "ddc/discrete_element.hpp"
#include "ddc/discrete_space.hpp"

namespace ddc {

template <class... DDim, std::enable_if_t<(sizeof...(DDim) > 1), int> = 0>
DDC_INLINE_FUNCTION Coordinate<typename DDim::continuous_dimension_type...> coordinate(
        DiscreteElement<DDim...> const& c)
{
    return Coordinate<typename DDim::continuous_dimension_type...>(coordinate(select<DDim>(c))...);
}

} // namespace ddc
