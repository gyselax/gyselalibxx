// SPDX-License-Identifier: MIT

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

    std::vector<double> m_derivs_vxmin_alloc;

    Span1D<double> m_derivs_vxmin;

    std::vector<double> m_derivs_vxmax_alloc;

    Span1D<double> m_derivs_vxmax;

public:
    SplineInterpolatorVx(
            SplineVxBuilder const& builder,
            SplineEvaluator<BSplinesVx> const& evaluator);

    ~SplineInterpolatorVx() override = default;

    DSpanVx operator()(DSpanVx inout_data, ViewVx<CoordVx> coordinates) const override;
};

class PreallocatableSplineInterpolatorVx : public IPreallocatableInterpolatorVx
{
    SplineVxBuilder const& m_builder;

    SplineEvaluator<BSplinesVx> const& m_evaluator;

public:
    PreallocatableSplineInterpolatorVx(
            SplineVxBuilder const& builder,
            SplineEvaluator<BSplinesVx> const& evaluator);

    ~PreallocatableSplineInterpolatorVx() override = default;

    InterpolatorVxProxy preallocate() const override;

    DSpanVx operator()(DSpanVx inout_data, ViewVx<CoordVx> coordinates) const override;
};
