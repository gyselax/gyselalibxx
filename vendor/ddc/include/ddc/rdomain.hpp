#pragma once

#include "ddc/rcoord.hpp"

template <class... RDim>
struct RDomain
{
    RCoord<RDim...> start;

    RCoord<RDim...> end;
};
