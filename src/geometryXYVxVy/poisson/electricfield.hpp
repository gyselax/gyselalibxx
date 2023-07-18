// SPDX-License-Identifier: MIT

#include <geometry.hpp>
#include <species_info.hpp>

//===========================================================================
// Compute efield = -dPhi/dx where Phi is the electrostatic potential
//  input : Phi values
//===========================================================================
class ElectricField
{
    SplineXYBuilder const& m_spline_xy_builder;

    SplineXYEvaluator m_spline_xy_evaluator;

public:
    ElectricField(
            SplineXYBuilder const& spline_xy_builder,
            SplineXYEvaluator const& spline_xy_evaluator);

    void operator()(
            DSpanXY electric_field_x,
            DSpanXY electric_field_y,
            DBSViewXY electrostatic_potential) const;
    void operator()(
            DSpanXY electric_field_x,
            DSpanXY electric_field_y,
            DViewXY electrostatic_potential) const;
};
