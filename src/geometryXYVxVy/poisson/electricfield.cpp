// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <sll/gauss_legendre_integration.hpp>

#include "electricfield.hpp"

ElectricField::ElectricField(
        SplineXYBuilder const& spline_xy_builder,
        SplineXYEvaluator const& spline_xy_evaluator)
    : m_spline_xy_builder(spline_xy_builder)
    , m_spline_xy_evaluator(spline_xy_evaluator)
{
}

//===========================================================================
// Compute efield = -dPhi/dx where Phi is the electrostatic potential
//  input : Phi values
//===========================================================================
void ElectricField::operator()(
        DSpanXY const electric_field_x,
        DSpanXY const electric_field_y,
        DBSViewXY const electrostatic_potential) const
{
    IDomainXY const& xy_dom = electric_field_x.domain();
    ddc::for_each(xy_dom, [&](IndexXY const ixy) {
        IndexX const ix = ddc::select<IDimX>(ixy);
        IndexY const iy = ddc::select<IDimY>(ixy);

        auto coordx = ddc::coordinate(ix);
        auto coordy = ddc::coordinate(iy);

        electric_field_x(ix, iy)
                = -m_spline_xy_evaluator
                           .deriv_dim_1(CoordXY(coordx, coordy), electrostatic_potential);
        electric_field_y(ix, iy)
                = -m_spline_xy_evaluator
                           .deriv_dim_2(CoordXY(coordx, coordy), electrostatic_potential);
    });
}


void ElectricField::operator()(
        DSpanXY const electric_field_x,
        DSpanXY const electric_field_y,
        DViewXY const electrostatic_potential) const
{
    Kokkos::Profiling::pushRegion("PoissonSolver");
    ddc::Chunk<double, BSDomainXY> elecpot_spline_coef((m_spline_xy_builder.spline_domain()));
    m_spline_xy_builder(elecpot_spline_coef, electrostatic_potential);
    (*this)(electric_field_x, electric_field_y, elecpot_spline_coef);
    Kokkos::Profiling::popRegion();
}
