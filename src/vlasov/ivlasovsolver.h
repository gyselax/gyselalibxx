#pragma once

class IVlasovSolver
{
public:
    virtual DBlockViewXVx& operator()(DBlockViewXVx& fdistribu, double mass_ratio, double dt) const = 0;
};
