#pragma once

class BoundaryValue
{
public:
    virtual double operator()(double x) const = 0;
};
