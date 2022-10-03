// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

namespace rdc {

struct DimX;

using DDimX = UniformPointSampling<DimX>;
using DElemX = DiscreteElement<DDimX>;
using DVectX = DiscreteVector<DDimX>;
using DDomX = DiscreteDomain<DDimX>;

void initialize_ddimx(Coordinate<DimX> origin, Coordinate<DimX> step);

} // namespace rdc
