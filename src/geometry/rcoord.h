#pragma once

#include "dim.h"
#include "taggedtuple.h"

using RCoordElement = double;

template <class... Tags>
using RCoord = TaggedTuple<RCoordElement, Tags...>;

using RCoordX = RCoord<Dim::X>;

using RCoordVx = RCoord<Dim::VX>;

using RCoordXVx = RCoord<Dim::X, Dim::VX>;
