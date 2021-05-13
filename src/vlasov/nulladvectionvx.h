#pragma once

#include "iadvectionvx.h"

class NullAdvectionVx : public IAdvectionVx
{
public:
    DBlockViewXVx operator()(DBlockViewXVx fdistribu, double mass_ratio, double dt) const override;
};
