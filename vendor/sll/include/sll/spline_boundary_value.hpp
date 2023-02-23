#pragma once
#include <functional>

#include <ddc/ddc.hpp>

#include <sll/polar_spline.hpp>

template <class BSplines>
class SplineBoundaryValue
{
public:
    virtual ~SplineBoundaryValue() = default;

    virtual double operator()(
            ddc::Coordinate<typename BSplines::tag_type> x,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplines>>) const = 0;
};

template <class BSplines1, class BSplines2>
class SplineBoundaryValue2D
{
public:
    virtual ~SplineBoundaryValue2D() = default;

    virtual double operator()(
            ddc::Coordinate<typename BSplines1::tag_type> x,
            ddc::Coordinate<typename BSplines2::tag_type> y,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplines1, BSplines2>>) const = 0;
};

template <class BSplines>
class PolarSplineBoundaryValue2D
{
public:
    virtual ~PolarSplineBoundaryValue2D() = default;

    virtual double operator()(double x, double y, PolarSplineView<BSplines>) const = 0;
};
