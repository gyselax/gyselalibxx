#pragma once

#include "rcoord.h"

template <class... Tags>
struct RDomain
{
    RCoord<Tags...> start;

    RCoord<Tags...> end;
};

using RDomainX = RDomain<Dim::X>;

using RDomainVx = RDomain<Dim::Vx>;

using RDomainXVx = RDomain<Dim::X, Dim::Vx>;
