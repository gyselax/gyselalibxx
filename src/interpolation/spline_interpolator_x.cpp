// SPDX-License-Identifier: MIT

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include "i_interpolator_x.hpp"
#include "spline_interpolator_x.hpp"

SplineInterpolatorX::SplineInterpolatorX(
        SplineXBuilder const& builder,
        SplineEvaluator<BSplinesX> const& evaluator)
    : m_builder(builder)
    , m_evaluator(evaluator)
    , m_coefs(builder.spline_domain())
{
}

DSpanX SplineInterpolatorX::operator()(DSpanX const inout_data, DViewX const coordinates) const
{
    m_builder(m_coefs, inout_data);
    m_evaluator(inout_data, coordinates, m_coefs);
    return inout_data;
}

PreallocatableSplineInterpolatorX::PreallocatableSplineInterpolatorX(
        SplineXBuilder const& builder,
        SplineEvaluator<BSplinesX> const& evaluator)
    : m_builder(builder)
    , m_evaluator(evaluator)
{
}

InterpolatorXProxy PreallocatableSplineInterpolatorX::preallocate() const
{
    return InterpolatorXProxy(std::make_unique<SplineInterpolatorX>(m_builder, m_evaluator));
}

DSpanX PreallocatableSplineInterpolatorX::operator()(
        DSpanX const inout_data,
        DViewX const coordinates) const
{
    return preallocate()(inout_data, coordinates);
}
