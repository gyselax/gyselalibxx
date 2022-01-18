#include <iostream>

#include <sll/gauss_legendre_integration.hpp>
#include <sll/spline_evaluator.hpp>

#include "electricfield.hpp"

ElectricField::ElectricField(
        SplineXBuilder const& spline_x_builder,
        SplineEvaluator<BSplinesX> const& spline_x_evaluator)
    : m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
{
}

//===========================================================================
// Compute efield = -dPhi/dx where Phi is the electrostatic potential
//  input : Phi values
//===========================================================================
void ElectricField::operator()(DSpanX electric_field, DBSViewX const electrostatic_potential) const
{
    IDomainX const& x_dom = electric_field.domain();
    for (IndexX const ix : x_dom) {
        electric_field(ix) = -m_spline_x_evaluator.deriv(to_real(ix), electrostatic_potential);
    }
}
void ElectricField::operator()(DSpanX electric_field, DViewX const electrostatic_potential) const
{
    Chunk<double, BSDomainX> elecpot_spline_coef(m_spline_x_builder.spline_domain());
    m_spline_x_builder(elecpot_spline_coef, electrostatic_potential);
    (*this)(electric_field, elecpot_spline_coef);
}
