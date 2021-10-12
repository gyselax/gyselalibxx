#pragma once

#include <vector>

#include <sll/spline_evaluator.h>

#include <geometry.h>

#include "iadvectionvx.h"

class BoundaryValue;
class SpeciesInformation;

class SplineAdvectionVx : public IAdvectionVx
{
private:
    const SplineVxBuilder& m_spline_vx_builder;

    SplineEvaluator<BSplinesVx> m_spline_vx_evaluator;

    std::vector<double> m_derivs_vxmin_data;

    Span1D<double> m_derivs_vxmin;

    std::vector<double> m_derivs_vxmax_data;

    Span1D<double> m_derivs_vxmax;

    SpeciesInformation const& m_species_info;

public:
    SplineAdvectionVx(
            SpeciesInformation const& species_info,
            const BSplinesVx& bspl,
            const SplineVxBuilder& spl_interp);

    SplineAdvectionVx(
            SpeciesInformation const& species_info,
            const BSplinesVx& bspl,
            const SplineVxBuilder& spl_interp,
            const BoundaryValue& bc_left,
            const BoundaryValue& bc_right);

    DSpanSpXVx operator()(DSpanSpXVx fdistribu, DViewX efield, double dt) const override;
};
