#pragma once

#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include "i_interpolator_vx.hpp"


class SplineInterpolatorVx : public IInterpolatorVx
{
private:
    SplineVxBuilder const& m_builder;

    SplineEvaluator<BSplinesVx> const& m_evaluator;

    mutable Chunk<double, BSDomainVx> m_coefs;

    std::vector<double> m_derivs_xmin_alloc;

    Span1D<double> m_derivs_xmin;

    std::vector<double> m_derivs_xmax_alloc;

    Span1D<double> m_derivs_xmax;

public:
    SplineInterpolatorVx(
            SplineVxBuilder const& builder,
            SplineEvaluator<BSplinesVx> const& evaluator);

    void operator()(DSpanVx const inout_data, DViewVx const coordinates) const override;
};

class PreallocatableSplineInterpolatorVx : public IPreallocatableInterpolatorVx
{
    SplineVxBuilder const& m_builder;

    SplineEvaluator<BSplinesVx> const& m_evaluator;

public:
    PreallocatableSplineInterpolatorVx(
            SplineVxBuilder const& builder,
            SplineEvaluator<BSplinesVx> const& evaluator);

    InterpolatorVxProxy preallocate() const override;

    void operator()(DSpanVx const inout_data, DViewVx const coordinates) const override;
};
