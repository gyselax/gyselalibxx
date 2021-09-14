#pragma once
#include "sll/boundary_value.h"

class NullBoundaryValue : public BoundaryValue
{
    NullBoundaryValue() = default;
    NullBoundaryValue(const NullBoundaryValue&) = delete;
    NullBoundaryValue(NullBoundaryValue&&) = delete;
    void operator=(const NullBoundaryValue&) = delete;
    void operator=(NullBoundaryValue&&) = delete;

public:
    inline virtual double operator()(double x) const override
    {
        return 0.0;
    }
    static NullBoundaryValue value;
};
