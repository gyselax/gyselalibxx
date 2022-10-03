#pragma once

#include "sll/spline_boundary_value.hpp"

template <class BSplines>
class NullBoundaryValue : public SplineBoundaryValue<BSplines>
{
public:
    NullBoundaryValue() = default;

    ~NullBoundaryValue() override = default;

    double operator()(double, ChunkSpan<const double, DiscreteDomain<BSplines>>) const final
    {
        return 0.0;
    }
};

template <class BSplines>
inline NullBoundaryValue<BSplines> const g_null_boundary;
