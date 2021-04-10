#pragma once

class IVlasovSolver
{
public:
    virtual void operator()(DBlock2D& cur, double mass_ratio, double dt) const = 0;
};
