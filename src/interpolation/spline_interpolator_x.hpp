#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include "i_interpolator_x.hpp"

class SplineInterpolatorX : public IInterpolatorX
{
private:
    SplineXBuilder const& m_builder;

    SplineEvaluator<BSplinesX> const& m_evaluator;

    mutable Chunk<double, BSDomainX> m_coefs;

public:
    SplineInterpolatorX(SplineXBuilder const& builder, SplineEvaluator<BSplinesX> const& evaluator);

    void operator()(DSpanX const inout_data, DViewX const coordinates) const override;
};

class PreallocatableSplineInterpolatorX : public IPreallocatableInterpolatorX
{
    SplineXBuilder const& m_builder;

    SplineEvaluator<BSplinesX> const& m_evaluator;

public:
    PreallocatableSplineInterpolatorX(
            SplineXBuilder const& builder,
            SplineEvaluator<BSplinesX> const& evaluator);

    InterpolatorXProxy preallocate() const override;

    void operator()(DSpanX const inout_data, DViewX const coordinates) const override;
};
