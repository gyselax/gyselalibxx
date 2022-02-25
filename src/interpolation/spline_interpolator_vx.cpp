// SPDX-License-Identifier: MIT

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include "i_interpolator_vx.hpp"
#include "spline_interpolator_vx.hpp"

SplineInterpolatorVx::SplineInterpolatorVx(
        SplineVxBuilder const& builder,
        SplineEvaluator<BSplinesVx> const& evaluator)
    : m_builder(builder)
    , m_evaluator(evaluator)
    , m_coefs(builder.spline_domain())
    , m_derivs_vxmin_alloc(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmin(m_derivs_vxmin_alloc.data(), m_derivs_vxmin_alloc.size())
    , m_derivs_vxmax_alloc(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmax(m_derivs_vxmax_alloc.data(), m_derivs_vxmax_alloc.size())
{
}

DSpanVx SplineInterpolatorVx::operator()(DSpanVx const inout_data, DViewVx const coordinates) const
{
    m_builder(m_coefs, inout_data, &m_derivs_vxmin, &m_derivs_vxmax);
    m_evaluator(inout_data, coordinates, m_coefs);
    return inout_data;
}

PreallocatableSplineInterpolatorVx::PreallocatableSplineInterpolatorVx(
        SplineVxBuilder const& builder,
        SplineEvaluator<BSplinesVx> const& evaluator)
    : m_builder(builder)
    , m_evaluator(evaluator)
{
}

InterpolatorVxProxy PreallocatableSplineInterpolatorVx::preallocate() const
{
    return InterpolatorVxProxy(std::make_unique<SplineInterpolatorVx>(m_builder, m_evaluator));
}

DSpanVx PreallocatableSplineInterpolatorVx::operator()(
        DSpanVx const inout_data,
        DViewVx const coordinates) const
{
    return preallocate()(inout_data, coordinates);
}
