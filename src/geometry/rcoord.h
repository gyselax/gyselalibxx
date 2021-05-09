#pragma once

#include "dim.h"
#include "taggedarray.h"

using RCoordElement = double;

template <class... Tags>
using RCoord = TaggedArray<RCoordElement, Tags...>;

using RCoordT = RCoord<Dim::T>;

using RCoordX = RCoord<Dim::X>;

using RCoordVx = RCoord<Dim::Vx>;

using RCoordXVx = RCoord<Dim::X, Dim::Vx>;



using RLengthElement = double;

template <class... Tags>
using RLength = TaggedArray<RLengthElement, Tags...>;

using RLengthT = RLength<Dim::T>;

using RLengthX = RLength<Dim::X>;

using RLengthVx = RLength<Dim::Vx>;

using RLengthXVx = RLength<Dim::X, Dim::Vx>;
