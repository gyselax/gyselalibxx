#pragma once

#include "taggedarray.h"

using RCoordElement = double;

template <class... Tags>
using RCoord = TaggedArray<RCoordElement, Tags...>;

using RLengthElement = double;

template <class... Tags>
using RLength = TaggedArray<RLengthElement, Tags...>;
