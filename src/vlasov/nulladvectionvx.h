#pragma once

#include "geometry.h"
#include "iadvectionvx.h"

class NullAdvectionVx : public IAdvectionVx
{
public:
    DSpanXVx operator()(DSpanXVx fdistribu, double mass_ratio, double dt) const override;
};
