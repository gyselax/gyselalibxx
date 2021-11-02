#pragma once

#include <geometry.hpp>

#include "iadvectionvx.hpp"

class NullAdvectionVx : public IAdvectionVx
{
public:
    DSpanSpXVx operator()(
            DSpanSpXVx const allfdistribu,
            DViewX const electric_potential,
            double const dt) const override;
};
