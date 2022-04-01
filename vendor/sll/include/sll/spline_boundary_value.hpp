#pragma once
#include <functional>

#include <ddc/ddc.hpp>

template <class BSplines>
class SplineBoundaryValue
{
public:
    virtual ~SplineBoundaryValue() = default;

    virtual double operator()(double x, ChunkSpan<double const, DiscreteDomain<BSplines>>)
            const = 0;
};
