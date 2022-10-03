#pragma once

#include "sll/spline_boundary_value.hpp"

template <class BSplines>
class NullBoundaryValue : public SplineBoundaryValue<BSplines>
{
    NullBoundaryValue() = default;

public:
    ~NullBoundaryValue() override = default;

    double operator()(double, ChunkSpan<const double, DiscreteDomain<BSplines>>) const final
    {
        return 0.0;
    }

    static inline NullBoundaryValue<BSplines> value;
};
