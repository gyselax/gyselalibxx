#pragma once

#include <geometry.h>

#include "iadvectionvx.h"

class NullAdvectionVx : public IAdvectionVx
{
public:
    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_potential, double dt)
            const override;
};
