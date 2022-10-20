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

template <class BSplines1, class BSplines2>
class SplineBoundaryValue2D
{
public:
    virtual ~SplineBoundaryValue2D() = default;

    virtual double operator()(
            double x,
            double y,
            ChunkSpan<double const, DiscreteDomain<BSplines1, BSplines2>>) const = 0;
};
