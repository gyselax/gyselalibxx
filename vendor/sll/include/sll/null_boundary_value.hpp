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

template <class BSplines1, class BSplines2>
class NullBoundaryValue2D : public SplineBoundaryValue2D<BSplines1, BSplines2>
{
public:
    NullBoundaryValue2D() = default;

    ~NullBoundaryValue2D() override = default;

    double operator()(
            double x,
            double y,
            ChunkSpan<double const, DiscreteDomain<BSplines1, BSplines2>>) const final
    {
        return 0.0;
    }
};

template <class BSplines1, class BSplines2>
inline NullBoundaryValue2D<BSplines1, BSplines2> const g_null_boundary_2d;
