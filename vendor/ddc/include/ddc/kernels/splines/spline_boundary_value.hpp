#pragma once
#include <functional>

#include <ddc/ddc.hpp>

namespace ddc {
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

} // namespace ddc
