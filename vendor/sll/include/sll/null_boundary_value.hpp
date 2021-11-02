#pragma once
#include "sll/boundary_value.hpp"

class NullBoundaryValue : public BoundaryValue
{
    NullBoundaryValue() = default;
    NullBoundaryValue(NullBoundaryValue const&) = delete;
    NullBoundaryValue(NullBoundaryValue&&) = delete;
    void operator=(NullBoundaryValue const&) = delete;
    void operator=(NullBoundaryValue&&) = delete;

public:
    inline virtual double operator()(double x) const override
    {
        return 0.0;
    }
    static NullBoundaryValue value;
};
