#pragma once

#include <ddc/BlockSpan>

class IVlasovSolver
{
public:
    virtual DSpanXVx operator()(
            DSpanXVx fdistribu,
            DViewX efield,
            double sqrt_me_on_mspecies,
            double dt) const = 0;
};
