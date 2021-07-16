#pragma once

#include "taggedarray.h"

using RCoordElement = double;

template <class... RDim>
using RCoord = TaggedArray<RCoordElement, RDim...>;

using RLengthElement = double;

template <class... RDim>
using RLength = TaggedArray<RLengthElement, RDim...>;
