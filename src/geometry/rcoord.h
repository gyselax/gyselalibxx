#pragma once

#include "dim.h"
#include "taggedarray.h"

using RCoordElement = double;

template <class... Tags>
using RCoord = TaggedArray<RCoordElement, Tags...>;

using RCoordX = RCoord<Dim::X>;

using RCoordVx = RCoord<Dim::Vx>;

using RCoordXVx = RCoord<Dim::X, Dim::Vx>;

using RLengthElement = double;

template <class... Tags>
using RLength = TaggedArray<RLengthElement, Tags...>;

using RLengthX = RLength<Dim::X>;

using RLengthVx = RLength<Dim::Vx>;

using RLengthXVx = RLength<Dim::X, Dim::Vx>;
