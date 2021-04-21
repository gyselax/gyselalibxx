#pragma once

#include <cstddef>

#include "dim.h"
#include "taggedarray.h"

using MCoordElement = std::size_t;

template <class... Tags>
using MCoord = TaggedArray<MCoordElement, Tags...>;

using MCoordX = MCoord<Dim::X>;

using MCoordVx = MCoord<Dim::Vx>;

using MCoordXVx = MCoord<Dim::X, Dim::Vx>;
