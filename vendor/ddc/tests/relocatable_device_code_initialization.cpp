// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include "relocatable_device_code_initialization.hpp"

namespace rdc {

void initialize_ddimx(Coordinate<DimX> const origin, Coordinate<DimX> const step)
{
    init_discrete_space<DDimX>(origin, step);
}

} // namespace rdc
