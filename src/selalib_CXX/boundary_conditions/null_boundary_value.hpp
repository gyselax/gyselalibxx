#pragma once
#include "boundary_value.hpp"

class NullBoundaryValue: public BoundaryValue
{
    public:
        NullBoundaryValue() = default;
        inline virtual double operator()(double x) const override {return 0.0; }
};