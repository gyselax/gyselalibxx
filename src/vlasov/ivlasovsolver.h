#pragma once

#include "block_span.h"

class IVlasovSolver
{
public:
    virtual DSpanXVx operator()(DSpanXVx fdistribu, double mass_ratio, double dt) const = 0;
};
