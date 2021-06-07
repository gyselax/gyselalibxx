#pragma once

#include "rcoord.h"

template <class... Tags>
struct RDomain
{
    RCoord<Tags...> start;

    RCoord<Tags...> end;
};
