#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>
#include <species_info.hpp>


//===========================================================================
// Compute efield = -dPhi/dx where Phi is the electrostatic potential
//  input : Phi values
//===========================================================================
class ElectricField
{
    SplineXBuilder const& m_spline_x_builder;

    SplineEvaluator<BSplinesX> m_spline_x_evaluator;

public:
    ElectricField(
            SplineXBuilder const& spline_x_builder,
            SplineEvaluator<BSplinesX> const& spline_x_evaluator);

    DSpanX operator()(DSpanX electric_field, DViewX const electrostatic_potential) const;

    DSpanX operator()(DSpanX electric_field, DBSViewX const electrostatic_potential) const;
};
