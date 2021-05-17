#pragma once

#include "iadvectionvx.h"

class NullAdvectionVx : public IAdvectionVx
{
public:
    DBlockSpanXVx operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt) const override;
};
