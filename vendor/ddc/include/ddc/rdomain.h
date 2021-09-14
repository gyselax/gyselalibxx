#pragma once

#include "ddc/rcoord.h"

template <class... RDim>
struct RDomain
{
    RCoord<RDim...> start;

    RCoord<RDim...> end;
};
