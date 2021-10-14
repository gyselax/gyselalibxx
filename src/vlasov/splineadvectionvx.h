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
    const SplineXBuilder& m_spline_x_builder;

    SplineEvaluator<BSplinesX> m_spline_x_evaluator;

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
            const BSplinesX& bspl_x,
            const SplineXBuilder& spl_x_interp,
            const BSplinesVx& bspl_vx,
            const SplineVxBuilder& spl_vx_interp);

    SplineAdvectionVx(
            SpeciesInformation const& species_info,
            const BSplinesX& bspl_x,
            const SplineXBuilder& spl_x_interp,
            const BSplinesVx& bspl_vx,
            const SplineVxBuilder& spl_vx_interp,
            const BoundaryValue& bc_vx_left,
            const BoundaryValue& bc_vx_right);

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_potential, double dt)
            const override;
};
