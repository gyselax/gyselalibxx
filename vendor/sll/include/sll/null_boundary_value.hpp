#pragma once
#include "sll/spline_boundary_value.hpp"

template <class BSplines>
class NullBoundaryValue : public SplineBoundaryValue<BSplines>
{
    NullBoundaryValue() = default;
    NullBoundaryValue(NullBoundaryValue const&) = delete;
    NullBoundaryValue(NullBoundaryValue&&) = delete;
    void operator=(NullBoundaryValue const&) = delete;
    void operator=(NullBoundaryValue&&) = delete;

public:
    ~NullBoundaryValue() override = default;

    inline double operator()(double x, ChunkSpan<const double, DiscreteDomain<BSplines>>)
            const final
    {
        return 0.0;
    }
    static inline NullBoundaryValue<BSplines> value;
};
