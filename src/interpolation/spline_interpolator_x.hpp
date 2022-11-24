// SPDX-License-Identifier: MIT

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

    ~SplineInterpolatorX() override = default;

    DSpanX operator()(DSpanX inout_data, ViewX<CoordX> coordinates) const override;
};

class PreallocatableSplineInterpolatorX : public IPreallocatableInterpolatorX
{
    SplineXBuilder const& m_builder;

    SplineEvaluator<BSplinesX> const& m_evaluator;

public:
    PreallocatableSplineInterpolatorX(
            SplineXBuilder const& builder,
            SplineEvaluator<BSplinesX> const& evaluator);

    ~PreallocatableSplineInterpolatorX() override = default;

    InterpolatorXProxy preallocate() const override;

    DSpanX operator()(DSpanX inout_data, ViewX<CoordX> coordinates) const override;
};
