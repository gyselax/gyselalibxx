#pragma once

#include <geometry.h>

#include "iadvectionvx.h"

class NullAdvectionVx : public IAdvectionVx
{
public:
    DSpanXVx operator()(
            DSpanXVx fdistribu,
            DViewX efield,
            int charge_species,
            double sqrt_me_on_mspecies,
            double dt) const override;
};
