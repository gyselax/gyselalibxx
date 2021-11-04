#pragma once

#include <geometry.hpp>

#include "iadvectionvx.hpp"

class NullAdvectionVx : public IAdvectionVx
{
public:
    NullAdvectionVx() = default;

    ~NullAdvectionVx() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_potential, double dt)
            const override;
};
