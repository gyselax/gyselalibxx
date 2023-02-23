#pragma once

#include "sll/spline_boundary_value.hpp"

template <class BSplines>
class NullBoundaryValue : public SplineBoundaryValue<BSplines>
{
public:
    NullBoundaryValue() = default;

    ~NullBoundaryValue() override = default;

    double operator()(
            ddc::Coordinate<typename BSplines::tag_type>,
            ddc::ChunkSpan<const double, ddc::DiscreteDomain<BSplines>>) const final
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
            ddc::Coordinate<typename BSplines1::tag_type> x,
            ddc::Coordinate<typename BSplines2::tag_type> y,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplines1, BSplines2>>) const final
    {
        return 0.0;
    }
};

template <class BSplines1, class BSplines2>
inline NullBoundaryValue2D<BSplines1, BSplines2> const g_null_boundary_2d;

template <class PolarBSplines>
class PolarNullBoundaryValue2D : public PolarSplineBoundaryValue2D<PolarBSplines>
{
public:
    PolarNullBoundaryValue2D() = default;

    ~PolarNullBoundaryValue2D() override = default;

    double operator()(double x, double y, PolarSplineView<PolarBSplines>) const final
    {
        return 0.0;
    }
};

template <class PolarBSplines>
inline PolarNullBoundaryValue2D<PolarBSplines> const g_polar_null_boundary_2d;
