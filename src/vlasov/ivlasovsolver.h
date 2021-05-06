#pragma once

class IVlasovSolver
{
public:
    virtual void operator()(DBlockViewXVx& cur, double mass_ratio, double dt) const = 0;
};
