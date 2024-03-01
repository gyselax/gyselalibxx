// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "electricfield.hpp"

ElectricField::ElectricField(
        SplineXBuilder_1d const& spline_x_builder,
        SplineXEvaluator_1d const& spline_x_evaluator)
    : m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
{
}

//===========================================================================
// Compute efield = -dPhi/dx where Phi is the electrostatic potential
//  input : Phi values
//===========================================================================
void ElectricField::operator()(
        host_t<DSpanX> const electric_field,
        host_t<DBSViewX> const electrostatic_potential) const
{
    IDomainX const& x_dom = electric_field.domain();
    ddc::for_each(x_dom, [&](IndexX const ix) {
        electric_field(ix)
                = -m_spline_x_evaluator.deriv(ddc::coordinate(ix), electrostatic_potential);
    });
}

void ElectricField::operator()(
        host_t<DSpanX> const electric_field,
        host_t<DViewX> const electrostatic_potential) const
{
    Kokkos::Profiling::pushRegion("ElectricField");
    host_t<ddc::Chunk<double, BSDomainX>> elecpot_spline_coef(m_spline_x_builder.bsplines_domain());
    m_spline_x_builder(elecpot_spline_coef.span_view(), electrostatic_potential);
    (*this)(electric_field, elecpot_spline_coef);
    Kokkos::Profiling::popRegion();
}
