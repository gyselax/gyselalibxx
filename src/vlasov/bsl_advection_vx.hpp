#pragma once

#include <vector>

#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>

#include "iadvectionvx.hpp"

class IPreallocatableInterpolatorX;
class IPreallocatableInterpolatorVx;
class BoundaryValue;
class SpeciesInformation;

class BslAdvectionVx : public IAdvectionVx
{
private:
    SplineXBuilder const& m_spline_x_builder;

    SplineEvaluator<BSplinesX> const& m_spline_x_evaluator;

    IPreallocatableInterpolatorVx const& m_interpolator_vx;

    SpeciesInformation const& m_species_info;

public:
    BslAdvectionVx(
            SpeciesInformation const& species_info,
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator,
            IPreallocatableInterpolatorVx const& interpolator_vx);

    ~BslAdvectionVx() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_potential, double dt)
            const override;
};
