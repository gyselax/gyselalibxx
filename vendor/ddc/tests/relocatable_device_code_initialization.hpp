// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

namespace rdc {

struct DimX;

using DDimX = ddc::UniformPointSampling<DimX>;
using DElemX = ddc::DiscreteElement<DDimX>;
using DVectX = ddc::DiscreteVector<DDimX>;
using DDomX = ddc::DiscreteDomain<DDimX>;

void initialize_ddimx(ddc::Coordinate<DimX> origin, ddc::Real step);

} // namespace rdc
