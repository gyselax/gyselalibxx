#pragma once

#include "geometry.h"
#include "iadvectionvx.h"

class NullAdvectionVx : public IAdvectionVx
{
public:
    DBlockSpanXVx operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt) const override;
};
