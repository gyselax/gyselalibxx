#pragma once

#include <array>
#include <cstddef>

#include "taggedtuple.h"

using MCoordElement = ptrdiff_t;

template <int NDIMS>
using MCoordND = std::array<MCoordElement, NDIMS>;

using MCoord1D = MCoordND<1>;

using MCoord2D = MCoordND<2>;

using MCoord3D = MCoordND<3>;


