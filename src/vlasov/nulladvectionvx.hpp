#pragma once

#include <geometry.hpp>

#include "iadvectionvx.hpp"

class NullAdvectionVx : public IAdvectionVx
{
public:
    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_potential, double dt)
            const override;
};
