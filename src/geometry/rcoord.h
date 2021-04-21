#pragma once

#include "dim.h"
#include "taggedarray.h"

using RCoordElement = double;

template <class... Tags>
using RCoord = TaggedArray<RCoordElement, Tags...>;

using RCoordX = RCoord<Dim::X>;

using RCoordVx = RCoord<Dim::Vx>;

using RCoordXVx = RCoord<Dim::X, Dim::Vx>;
